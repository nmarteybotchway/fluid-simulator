//
// Created by neil on 8/30/25.
//

#include "Fluid.h"

Fluid::Fluid(int size, int iters, float dt, float diffusion, float viscosity) :
    N{size}, iters{iters}, dt{dt}, diff{diffusion}, visc{viscosity} {
    s.resize(N*N, 0.0f);
    density.resize(N*N, 0.0f);

    Vx.resize(N*N, 0.0f);
    Vy.resize(N*N, 0.0f);

    Vx0.resize(N*N, 0.0f);
    Vy0.resize(N*N, 0.0f);
}

inline int Fluid::IX(int x, int y) const {
    return x + y * N;
}

void Fluid::addDensity(int x, int y, float amount) {
    density[IX(x, y)] += amount;
}

void Fluid::addVelocity(int x, int y, float amountX, float amountY) {
    Vx[IX(x, y)] += amountX;
    Vy[IX(x, y)] += amountY;
}

void Fluid::set_bnd(int b, std::vector<float>& x) {
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

void Fluid::lin_solve(int b, std::vector<float>& x, const std::vector<float>& x0, float a, float c) {
    float cRecip = 1.0f / c;
    for (int k = 0; k < iters; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                x[IX(i, j)] =
                    (x0[IX(i, j)]
                    + a*(x[IX(i+1, j)] + x[IX(i-1, j)]
                         + x[IX(i, j+1)] + x[IX(i, j-1)])) * cRecip;
            }
        }
        set_bnd(b, x);
    }
}

void Fluid::project(std::vector<float>& velocX, std::vector<float>& velocY,
             std::vector<float>& p, std::vector<float>& div) {
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            div[IX(i, j)] = -0.5f*(
                velocX[IX(i+1, j)] - velocX[IX(i-1, j)]
                + velocY[IX(i, j+1)] - velocY[IX(i, j-1)]
            ) / N;
            p[IX(i, j)] = 0;
        }
    }
    set_bnd(0, div);
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 6);

    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            velocX[IX(i, j)] -= 0.5f * (p[IX(i+1, j)] - p[IX(i-1, j)]) * N;
            velocY[IX(i, j)] -= 0.5f * (p[IX(i, j+1)] - p[IX(i, j-1)]) * N;
        }
    }
    set_bnd(1, velocX);
    set_bnd(2, velocY);
}

void Fluid::step() {
    diffuse(1, Vx0, Vx, visc, dt);
    diffuse(2, Vy0, Vy, visc, dt);

    project(Vx0, Vy0, Vx, Vy);

    advect(1, Vx, Vx0, Vx0, Vy0, dt);
    advect(2, Vy, Vy0, Vx0, Vy0, dt);

    project(Vx, Vy, Vx0, Vy0);

    diffuse(0, s, density, diff, dt);
    advect(0, density, s, Vx, Vy, dt);
}

void Fluid::advect(int b, std::vector<float>& d, const std::vector<float>& d0,
            const std::vector<float>& velocX, const std::vector<float>& velocY, float dt) {
    float Nfloat = N;
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            float x = i - dtx * velocX[IX(i,j)];
            float y = j - dty * velocY[IX(i,j)];

            if (x < 0.5f) x = 0.5f;
            if (x > Nfloat - 1.5f) x = Nfloat - 1.5f;
            if (y < 0.5f) y = 0.5f;
            if (y > Nfloat - 1.5f) y = Nfloat - 1.5f;

            int i0 = static_cast<int>(floor(x));
            int i1 = i0 + 1;
            int j0 = static_cast<int>(floor(y));
            int j1 = j0 + 1;

            float s1 = x - i0;
            float s0 = 1.0f - s1;
            float t1 = y - j0;
            float t0 = 1.0f - t1;

            d[IX(i,j)] =
                s0*(t0*d0[IX(i0,j0)] + t1*d0[IX(i0,j1)]) +
                s1*(t0*d0[IX(i1,j0)] + t1*d0[IX(i1,j1)]);
        }
    }
    set_bnd(b, d);
}

void Fluid::diffuse(int b, std::vector<float>& x, const std::vector<float>& x0,
             float diff, float dt) {
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6*a);
}

void Fluid::renderD(sf::RenderWindow& window) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            float d = density[IX(i,j)];
            sf::Uint8 c = static_cast<sf::Uint8>(std::clamp(d, 0.0f, 255.0f));
            sf::RectangleShape cell(sf::Vector2f(10,10));
            cell.setPosition(i*10, j*10);
            cell.setFillColor(sf::Color(c,c,c));
            window.draw(cell);
        }
    }
}

void Fluid::fadeD() {
    for (int i = 0; i < this->density.size(); i++) {
        float d = density[i];
        density[i] = std::clamp(d-0.05f, 0.f, 255.f);
    }
}