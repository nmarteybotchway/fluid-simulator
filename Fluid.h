#ifndef FLUID_SIMULATION_FLUID_H
#define FLUID_SIMULATION_FLUID_H

#include <cmath>
#include <vector>
#include <SFML/Graphics.hpp>
#include <algorithm>

class Fluid {
    int gridSize;
    int solverIterations;
    float timeStep;
    float diffusionRate;
    float viscosity;
    int scale;
    std::vector<float> prevDensity;
    std::vector<float> density;

    std::vector<float> velocityX;
    std::vector<float> velocityY;

    std::vector<float> velocityXPrev;
    std::vector<float> velocityYPrev;

public:
    Fluid(int size, int iters, float dt, float diffusion, float viscosity, int scale);
    [[nodiscard]] int getScale() const { return scale; }
    [[nodiscard]] int getGridSize() const { return gridSize; }
    [[nodiscard]] inline int IX(int x, int y) const;
    [[nodiscard]] int getWindowWidth() const { return gridSize * scale; }
    [[nodiscard]] int getWindowHeight() const { return gridSize * scale; }
    void addDensity(int x, int y, float amount);
    void addVelocity(int x, int y, float amountX, float amountY);

    void applyBoundaryConditions(int b, std::vector<float>& x) const;
    void linearSolve(int b, std::vector<float>& x, const std::vector<float>& x0, float a, float c) const;
    void projectVelocity(std::vector<float>& velocityX, std::vector<float>& velocityY,
                 std::vector<float>& p, std::vector<float>& div) const;
    void advanceSimulation();

    void advectQuantity(int b, std::vector<float>& d, const std::vector<float>& d0,
                const std::vector<float>& velocityX, const std::vector<float>& velocityY, float dt) const;

    void diffuseQuantity(int b, std::vector<float>& x, const std::vector<float>& x0,
                 float diff, float dt) const;
    void renderDensity(sf::RenderWindow& window) const;
    void decayDensity();
};

#endif //FLUID_SIMULATION_FLUID_H
