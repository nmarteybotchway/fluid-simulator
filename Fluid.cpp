#include "Fluid.h"
#include <cmath>

Fluid::Fluid(const int gridSize, const int solverIterations, const float timeStep, const float diffusionRate,
    const float viscosity, const int scale) :
    gridSize{gridSize},
    solverIterations{solverIterations},
    timeStep{timeStep},
    diffusionRate{diffusionRate},
    viscosity{viscosity},
    scale{scale},
    prevDensity{gridSize},
    density{gridSize},
    velocityX{gridSize},
    velocityY{gridSize},
    velocityXPrev{gridSize},
    velocityYPrev{gridSize} {}

void Fluid::addDensity(const int x, const int y, const float amount) {
    constexpr int radius = 1;  // adjust this for smoother / wider strokes

    for (int j = -radius; j <= radius; j++) {
        for (int i = -radius; i <= radius; i++) {
            const int xi = x + i;
            int yj = y + j;
            if (xi >= 0 && xi < gridSize && yj >= 0 && yj < gridSize) {
                const float dist2 = i*i + j*j;
                const float weight = std::exp(-dist2 / (2.0f * radius * radius));
                density(xi, yj) += amount * weight;
            }
        }
    }
}

void Fluid::addVelocity(const int x, const int y, const float amountX, const float amountY) {
    velocityX(x, y) += amountX;
    velocityY(x, y) += amountY;
}

void Fluid::applyBoundaryConditions(const int b, Grid& x) const {
    for (int i = 1; i < gridSize - 1; i++) {
        x(i, 0) = b == 2 ? -x(i, 1) : x(i, 1);
        x(i, gridSize - 1) = b == 2 ? -x(i, gridSize - 2) : x(i, gridSize - 2);
    }
    for (int j = 1; j < gridSize - 1; j++) {
        x(0, j) = b == 1 ? -x(1, j) : x(1, j);
        x(gridSize - 1, j) = b == 1 ? -x(gridSize - 2, j) : x(gridSize - 2, j);
    }

    x(0, 0) = 0.5f * (x(1, 0) + x(0, 1));
    x(0, gridSize-1) = 0.5f * (x(1, gridSize-1) + x(0, gridSize-2));
    x(gridSize-1, 0) = 0.5f * (x(gridSize-2, 0) + x(gridSize-1, 1));
    x(gridSize-1, gridSize-1) = 0.5f * (x(gridSize-2, gridSize-1) + x(gridSize-1, gridSize-2));
}

void Fluid::linearSolve(const int b, Grid& x, const Grid& x0, const float a, const float c)
const {
    const float cRecip = 1.0f / c;
    Grid xNew(x.size, 0.0f);

    for (int k = 0; k < solverIterations; k++) {
#pragma omp parallel for collapse(2)
        for (int j = 1; j < gridSize - 1; j++) {
            for (int i = 1; i < gridSize - 1; i++) {
                xNew(i, j) =
                    (x0(i, j)
                    + a*(x(i+1, j) + x(i-1, j)
                         + x(i, j+1) + x(i, j-1))) * cRecip;
            }
        }

        // swap roles: new becomes current
        std::swap(x, xNew);

        // apply boundary conditions
        applyBoundaryConditions(b, x);
    }
}

void Fluid::projectVelocity(Grid& velocityX, Grid& velocityY,
                    Grid& p, Grid& div) const {
    // First loop: compute divergence & reset p
#pragma omp parallel for collapse(2)
    for (int j = 1; j < gridSize - 1; j++) {
        for (int i = 1; i < gridSize - 1; i++) {
            div(i, j) = -0.5f*(
                velocityX(i+1, j) - velocityX(i-1, j)
                + velocityY(i, j+1) - velocityY(i, j-1)
            ) / gridSize;
            p(i, j) = 0;
        }
    }

    applyBoundaryConditions(0, div);
    applyBoundaryConditions(0, p);
    linearSolve(0, p, div, 1, 6);

    // Second loop: subtract gradient of pressure
#pragma omp parallel for collapse(2)
    for (int j = 1; j < gridSize - 1; j++) {
        for (int i = 1; i < gridSize - 1; i++) {
            velocityX(i, j) -= 0.5f * (p(i+1, j) - p(i-1, j)) * gridSize;
            velocityY(i, j) -= 0.5f * (p(i, j+1) - p(i, j-1)) * gridSize;
        }
    }

    applyBoundaryConditions(1, velocityX);
    applyBoundaryConditions(2, velocityY);
}


void Fluid::advanceSimulation() {
    // Diffuse velocities
#pragma omp parallel sections
    {
#pragma omp section
        diffuseQuantity(1, velocityXPrev, velocityX, viscosity, timeStep);

#pragma omp section
        diffuseQuantity(2, velocityYPrev, velocityY, viscosity, timeStep);
    }

    projectVelocity(velocityXPrev, velocityYPrev, velocityX, velocityY);

    // Advect velocities
#pragma omp parallel sections
    {
#pragma omp section
        advectQuantity(1, velocityX, velocityXPrev, velocityXPrev, velocityYPrev, timeStep);

#pragma omp section
        advectQuantity(2, velocityY, velocityYPrev, velocityXPrev, velocityYPrev, timeStep);
    }

    projectVelocity(velocityX, velocityY, velocityXPrev, velocityYPrev);

    // Diffuse and advect density
    diffuseQuantity(0, prevDensity, density, diffusionRate, timeStep);
    advectQuantity(0, density, prevDensity, velocityX, velocityY, timeStep);
}

void Fluid::advectQuantity(int b, Grid& d, const Grid& d0,
                   const Grid& velocityX, const Grid& velocityY, float timeStep) const {
    float nfloat = gridSize;
    float dtx = timeStep * (gridSize - 2);
    float dty = timeStep * (gridSize - 2);

    // Parallelize the 2D grid loop
#pragma omp parallel for collapse(2)
    for (int j = 1; j < gridSize - 1; j++) {
        for (int i = 1; i < gridSize - 1; i++) {
            float x = i - dtx * velocityX(i,j);
            float y = j - dty * velocityY(i,j);

            if (x < 0.5f) x = 0.5f;
            if (x > nfloat - 1.5f) x = nfloat - 1.5f;
            if (y < 0.5f) y = 0.5f;
            if (y > nfloat - 1.5f) y = nfloat - 1.5f;

            const int i0 = static_cast<int>(std::floor(x));
            const int i1 = i0 + 1;
            const int j0 = static_cast<int>(std::floor(y));
            const int j1 = j0 + 1;

            const float s1 = x - i0;
            const float s0 = 1.0f - s1;
            const float t1 = y - j0;
            const float t0 = 1.0f - t1;

            d(i,j) =
                s0*(t0*d0(i0,j0) + t1*d0(i0,j1)) +
                s1*(t0*d0(i1,j0) + t1*d0(i1,j1));
        }
    }

    applyBoundaryConditions(b, d);
}

void Fluid::diffuseQuantity(int b, Grid& x, const Grid& x0,
             const float diffusionRate, const float timeStep) const {
    const float a = timeStep * diffusionRate * (gridSize - 2) * (gridSize - 2);
    linearSolve(b, x, x0, a, 1 + 6*a);
}

void Fluid::renderDensity(sf::RenderWindow& window) const {
    if (!densityTexInit) {
        densityTex.create(gridSize, gridSize);
        densityTex.setSmooth(true);
        densityTexInit = true;
    }


    // Update texture pixels from density
    std::vector<sf::Uint8> pixels(gridSize * gridSize * 4); // RGBA
    for (int j = 0; j < gridSize; j++) {
        for (int i = 0; i < gridSize; i++) {
            const float d = density(i, j) / 255.0f;
            const sf::Uint8 c = static_cast<sf::Uint8>(std::clamp(d * 255.0f, 0.f, 255.f));
            const int idx = 4 * (j * gridSize + i);
            pixels[idx + 0] = c; // R
            pixels[idx + 1] = c; // G
            pixels[idx + 2] = c; // B
            pixels[idx + 3] = 255; // A
        }
    }

    densityTex.update(pixels.data());

    sf::Sprite sprite(densityTex);
    sprite.setScale(scale, scale); // use scale member

    window.draw(sprite);
}

void Fluid::decayDensity() {
    constexpr float densityDecayRate = 0.05f;

#pragma omp parallel for collapse(2)
    for (int y = 0; y < gridSize; y++) {
        for (int x = 0; x < gridSize; x++) {
            float& d = density(x, y);
            d = std::clamp(d - densityDecayRate, 0.0f, 255.0f);
        }
    }
}
