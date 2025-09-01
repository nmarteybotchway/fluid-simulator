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
    int scale;
    std::vector<float> s;
    std::vector<float> density;

    std::vector<float> Vx;
    std::vector<float> Vy;

    std::vector<float> Vx0;
    std::vector<float> Vy0;

public:
    Fluid(int size, int iters, float dt, float diffusion, float viscosity, int scale);
    int getScale() { return scale; }
    inline int IX(int x, int y) const;

    void addDensity(int x, int y, float amount);
    void addVelocity(int x, int y, float amountX, float amountY);

    void set_bnd(int b, std::vector<float>& x) const;
    void lin_solve(int b, std::vector<float>& x, const std::vector<float>& x0, float a, float c) const;
    void project(std::vector<float>& velocityX, std::vector<float>& velocityY,
                 std::vector<float>& p, std::vector<float>& div) const;
    void step();

    void advect(int b, std::vector<float>& d, const std::vector<float>& d0,
                const std::vector<float>& velocityX, const std::vector<float>& velocityY, float dt) const;

    void diffuse(int b, std::vector<float>& x, const std::vector<float>& x0,
                 float diff, float dt) const;
    void renderD(sf::RenderWindow& window) const;
    void fadeD();
};

#endif //FLUID_SIMULATION_FLUID_H
