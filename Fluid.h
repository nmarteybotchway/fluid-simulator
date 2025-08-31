#ifndef FLUID_SIMULATION_FLUID_H
#define FLUID_SIMULATION_FLUID_H

#include <cmath>
#include <vector>
#include <SFML/Graphics.hpp>
#include <algorithm>

class Fluid {
    int N;
    int iters;
    float dt;
    float diff;
    float visc;

    std::vector<float> s;
    std::vector<float> density;

    std::vector<float> Vx;
    std::vector<float> Vy;

    std::vector<float> Vx0;
    std::vector<float> Vy0;

public:
    Fluid(int size, int iters, float dt, float diffusion, float viscosity);

    inline int IX(int x, int y) const;

    void addDensity(int x, int y, float amount);
    void addVelocity(int x, int y, float amountX, float amountY);

    void set_bnd(int b, std::vector<float>& x);
    void lin_solve(int b, std::vector<float>& x, const std::vector<float>& x0, float a, float c);
    void project(std::vector<float>& velocX, std::vector<float>& velocY,
                 std::vector<float>& p, std::vector<float>& div);
    void step();

    void advect(int b, std::vector<float>& d, const std::vector<float>& d0,
                const std::vector<float>& velocX, const std::vector<float>& velocY, float dt);

    void diffuse(int b, std::vector<float>& x, const std::vector<float>& x0,
                 float diff, float dt);
    void renderD(sf::RenderWindow& window);
    void fadeD();
};

#endif //FLUID_SIMULATION_FLUID_H
