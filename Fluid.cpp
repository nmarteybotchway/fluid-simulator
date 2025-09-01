//
// Created by neil on 8/30/25.
//

#include "Fluid.h"

#include <cmath>

Fluid::Fluid(const int gridSize, const int solverIterations, const float timeStep, const float diffusionRate, const float viscosity, const int scale) :
    gridSize{gridSize}, solverIterations{solverIterations}, timeStep{timeStep}, diffusionRate{diffusionRate}, viscosity{viscosity}, scale{scale}{
    prevDensity.resize(gridSize*gridSize, 0.0f);
    density.resize(gridSize*gridSize, 0.0f);

    velocityX.resize(gridSize*gridSize, 0.0f);
    velocityY.resize(gridSize*gridSize, 0.0f);

    velocityXPrev.resize(gridSize*gridSize, 0.0f);
    velocityYPrev.resize(gridSize*gridSize, 0.0f);
}

inline int Fluid::IX(const int x, const int y) const {
    return x + y * gridSize;
}

void Fluid::addDensity(const int x, const int y, const float amount) {
    constexpr int radius = 1;  // adjust this for smoother / wider strokes

    for (int j = -radius; j <= radius; j++) {
        for (int i = -radius; i <= radius; i++) {
            const int xi = x + i;
            int yj = y + j;
            if (xi >= 0 && xi < gridSize && yj >= 0 && yj < gridSize) {
                const float dist2 = i*i + j*j;
                const float weight = std::exp(-dist2 / (2.0f * radius * radius));
                density[IX(xi, yj)] += amount * weight;
            }
        }
    }
}

void Fluid::addVelocity(const int x, const int y, const float amountX, const float amountY) {
    velocityX[IX(x, y)] += amountX;
    velocityY[IX(x, y)] += amountY;
}

void Fluid::applyBoundaryConditions(const int b, std::vector<float>& x) const {
    for (int i = 1; i < gridSize - 1; i++) {
        x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, gridSize - 1)] = b == 2 ? -x[IX(i, gridSize - 2)] : x[IX(i, gridSize - 2)];
    }
    for (int j = 1; j < gridSize - 1; j++) {
        x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
        x[IX(gridSize - 1, j)] = b == 1 ? -x[IX(gridSize - 2, j)] : x[IX(gridSize - 2, j)];
    }

    x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, gridSize-1)] = 0.5f * (x[IX(1, gridSize-1)] + x[IX(0, gridSize-2)]);
    x[IX(gridSize-1, 0)] = 0.5f * (x[IX(gridSize-2, 0)] + x[IX(gridSize-1, 1)]);
    x[IX(gridSize-1, gridSize-1)] = 0.5f * (x[IX(gridSize-2, gridSize-1)] + x[IX(gridSize-1, gridSize-2)]);
}

void Fluid::linearSolve(const int b, std::vector<float>& x, const std::vector<float>& x0, const float a, const float c)
const {
    const float cRecip = 1.0f / c;
    std::vector<float> xNew(x.size(), 0.0f);

    for (int k = 0; k < solverIterations; k++) {
#pragma omp parallel for collapse(2)
        for (int j = 1; j < gridSize - 1; j++) {
            for (int i = 1; i < gridSize - 1; i++) {
                xNew[IX(i, j)] =
                    (x0[IX(i, j)]
                    + a*(x[IX(i+1, j)] + x[IX(i-1, j)]
                         + x[IX(i, j+1)] + x[IX(i, j-1)])) * cRecip;
            }
        }

        // swap roles: new becomes current
        std::swap(x, xNew);

        // apply boundary conditions
        applyBoundaryConditions(b, x);
    }
}

void Fluid::projectVelocity(std::vector<float>& velocityX, std::vector<float>& velocityY,
                    std::vector<float>& p, std::vector<float>& div) const {
    // First loop: compute divergence & reset p
#pragma omp parallel for collapse(2)
    for (int j = 1; j < gridSize - 1; j++) {
        for (int i = 1; i < gridSize - 1; i++) {
            div[IX(i, j)] = -0.5f*(
                velocityX[IX(i+1, j)] - velocityX[IX(i-1, j)]
                + velocityY[IX(i, j+1)] - velocityY[IX(i, j-1)]
            ) / gridSize;
            p[IX(i, j)] = 0;
        }
    }

    applyBoundaryConditions(0, div);
    applyBoundaryConditions(0, p);
    linearSolve(0, p, div, 1, 6);

    // Second loop: subtract gradient of pressure
#pragma omp parallel for collapse(2)
    for (int j = 1; j < gridSize - 1; j++) {
        for (int i = 1; i < gridSize - 1; i++) {
            velocityX[IX(i, j)] -= 0.5f * (p[IX(i+1, j)] - p[IX(i-1, j)]) * gridSize;
            velocityY[IX(i, j)] -= 0.5f * (p[IX(i, j+1)] - p[IX(i, j-1)]) * gridSize;
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

void Fluid::advectQuantity(int b, std::vector<float>& d, const std::vector<float>& d0,
                   const std::vector<float>& velocityX, const std::vector<float>& velocityY, float timeStep) const {
    float nfloat = gridSize;
    float dtx = timeStep * (gridSize - 2);
    float dty = timeStep * (gridSize - 2);

    // Parallelize the 2D grid loop
#pragma omp parallel for collapse(2)
    for (int j = 1; j < gridSize - 1; j++) {
        for (int i = 1; i < gridSize - 1; i++) {
            float x = i - dtx * velocityX[IX(i,j)];
            float y = j - dty * velocityY[IX(i,j)];

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

            d[IX(i,j)] =
                s0*(t0*d0[IX(i0,j0)] + t1*d0[IX(i0,j1)]) +
                s1*(t0*d0[IX(i1,j0)] + t1*d0[IX(i1,j1)]);
        }
    }

    applyBoundaryConditions(b, d);
}

void Fluid::diffuseQuantity(int b, std::vector<float>& x, const std::vector<float>& x0,
             float diffusionRate, float timeStep) const {
    float a = timeStep * diffusionRate * (gridSize - 2) * (gridSize - 2);
    linearSolve(b, x, x0, a, 1 + 6*a);
}

void Fluid::renderDensity(sf::RenderWindow& window) const {
    static sf::Texture densityTex;
    static bool initialized = false;

    if (!initialized) {
        densityTex.create(gridSize, gridSize);
        densityTex.setSmooth(true); // smooth scaling
        initialized = true;
    }

    // Update texture pixels from density
    std::vector<sf::Uint8> pixels(gridSize * gridSize * 4); // RGBA
    for (int j = 0; j < gridSize; j++) {
        for (int i = 0; i < gridSize; i++) {
            const float d = density[IX(i, j)] / 255.0f;
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
#pragma omp parallel for
    for (int i = 0; i < this->density.size(); i++) {
        const float d = density[i];
        density[i] = std::clamp(d - 0.05f, 0.f, 255.f);
    }
}
