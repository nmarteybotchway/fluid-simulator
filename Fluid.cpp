//
// Created by neil on 8/30/25.
//

#include "Fluid.h"

#include <cmath>

Fluid::Fluid(const int size, const int iters, const float dt, const float diffusion, const float viscosity, const int scale) :
    N{size}, iters{iters}, dt{dt}, diff{diffusion}, visc{viscosity}, scale{scale}{
    s.resize(N*N, 0.0f);
    density.resize(N*N, 0.0f);

    Vx.resize(N*N, 0.0f);
    Vy.resize(N*N, 0.0f);

    Vx0.resize(N*N, 0.0f);
    Vy0.resize(N*N, 0.0f);
}

inline int Fluid::IX(const int x, const int y) const {
    return x + y * N;
}

void Fluid::addDensity(const int x, const int y, const float amount) {
    constexpr int radius = 1;  // adjust this for smoother / wider strokes

    for (int j = -radius; j <= radius; j++) {
        for (int i = -radius; i <= radius; i++) {
            const int xi = x + i;
            int yj = y + j;
            if (xi >= 0 && xi < N && yj >= 0 && yj < N) {
                const float dist2 = i*i + j*j;
                const float weight = std::exp(-dist2 / (2.0f * radius * radius));
                density[IX(xi, yj)] += amount * weight;
            }
        }
    }
}

void Fluid::addVelocity(const int x, const int y, const float amountX, const float amountY) {
    Vx[IX(x, y)] += amountX;
    Vy[IX(x, y)] += amountY;
}

void Fluid::set_bnd(const int b, std::vector<float>& x) const {
    for (int i = 1; i < N - 1; i++) {
        x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N - 1)] = b == 2 ? -x[IX(i, N - 2)] : x[IX(i, N - 2)];
    }
    for (int j = 1; j < N - 1; j++) {
        x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
        x[IX(N - 1, j)] = b == 1 ? -x[IX(N - 2, j)] : x[IX(N - 2, j)];
    }

    x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N-1)] = 0.5f * (x[IX(1, N-1)] + x[IX(0, N-2)]);
    x[IX(N-1, 0)] = 0.5f * (x[IX(N-2, 0)] + x[IX(N-1, 1)]);
    x[IX(N-1, N-1)] = 0.5f * (x[IX(N-2, N-1)] + x[IX(N-1, N-2)]);
}

void Fluid::lin_solve(const int b, std::vector<float>& x, const std::vector<float>& x0, const float a, const float c)
const {
    const float cRecip = 1.0f / c;
    std::vector<float> xNew(x.size(), 0.0f);

    for (int k = 0; k < iters; k++) {
#pragma omp parallel for collapse(2)
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                xNew[IX(i, j)] =
                    (x0[IX(i, j)]
                    + a*(x[IX(i+1, j)] + x[IX(i-1, j)]
                         + x[IX(i, j+1)] + x[IX(i, j-1)])) * cRecip;
            }
        }

        // swap roles: new becomes current
        std::swap(x, xNew);

        // apply boundary conditions
        set_bnd(b, x);
    }
}

void Fluid::project(std::vector<float>& velocityX, std::vector<float>& velocityY,
                    std::vector<float>& p, std::vector<float>& div) const {
    // First loop: compute divergence & reset p
#pragma omp parallel for collapse(2)
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            div[IX(i, j)] = -0.5f*(
                velocityX[IX(i+1, j)] - velocityX[IX(i-1, j)]
                + velocityY[IX(i, j+1)] - velocityY[IX(i, j-1)]
            ) / N;
            p[IX(i, j)] = 0;
        }
    }

    set_bnd(0, div);
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 6);

    // Second loop: subtract gradient of pressure
#pragma omp parallel for collapse(2)
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            velocityX[IX(i, j)] -= 0.5f * (p[IX(i+1, j)] - p[IX(i-1, j)]) * N;
            velocityY[IX(i, j)] -= 0.5f * (p[IX(i, j+1)] - p[IX(i, j-1)]) * N;
        }
    }

    set_bnd(1, velocityX);
    set_bnd(2, velocityY);
}


void Fluid::step() {
    // Diffuse velocities
#pragma omp parallel sections
    {
#pragma omp section
        diffuse(1, Vx0, Vx, visc, dt);

#pragma omp section
        diffuse(2, Vy0, Vy, visc, dt);
    }

    project(Vx0, Vy0, Vx, Vy);

    // Advect velocities
#pragma omp parallel sections
    {
#pragma omp section
        advect(1, Vx, Vx0, Vx0, Vy0, dt);

#pragma omp section
        advect(2, Vy, Vy0, Vx0, Vy0, dt);
    }

    project(Vx, Vy, Vx0, Vy0);

    // Diffuse and advect density
    diffuse(0, s, density, diff, dt);
    advect(0, density, s, Vx, Vy, dt);
}

void Fluid::advect(int b, std::vector<float>& d, const std::vector<float>& d0,
                   const std::vector<float>& velocityX, const std::vector<float>& velocityY, float dt) const {
    float nfloat = N;
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);

    // Parallelize the 2D grid loop
#pragma omp parallel for collapse(2)
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
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

    set_bnd(b, d);
}

void Fluid::diffuse(int b, std::vector<float>& x, const std::vector<float>& x0,
             float diff, float dt) const {
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6*a);
}

void Fluid::renderD(sf::RenderWindow& window) const {
    static sf::Texture densityTex;
    static bool initialized = false;

    if (!initialized) {
        densityTex.create(N, N);
        densityTex.setSmooth(true); // smooth scaling
        initialized = true;
    }

    // Update texture pixels from density
    std::vector<sf::Uint8> pixels(N * N * 4); // RGBA
    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            const float d = density[IX(i, j)] / 255.0f;
            const sf::Uint8 c = static_cast<sf::Uint8>(std::clamp(d * 255.0f, 0.f, 255.f));
            const int idx = 4 * (j * N + i);
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

void Fluid::fadeD() {
#pragma omp parallel for
    for (int i = 0; i < this->density.size(); i++) {
        const float d = density[i];
        density[i] = std::clamp(d - 0.05f, 0.f, 255.f);
    }
}
