# pragma once
#include <cmath>
#include <SFML/Graphics.hpp>
#include "Grid.h"

class Fluid {
    int gridSize;
    int solverIterations;
    float timeStep;
    float diffusionRate;
    float viscosity;
    int scale;
    Grid prevDensity;
    Grid density;
    Grid velocityX;
    Grid velocityY;
    Grid velocityXPrev;
    Grid velocityYPrev;
    mutable sf::Texture densityTex;
    mutable bool densityTexInit = false;

public:
    Fluid(int size, int iters, float dt, float diffusion, float viscosity, int scale);
    [[nodiscard]] int getScale() const { return scale; }
    [[nodiscard]] int getGridSize() const { return gridSize; }
    [[nodiscard]] int getWindowWidth() const { return gridSize * scale; }
    [[nodiscard]] int getWindowHeight() const { return gridSize * scale; }
    void addDensity(int x, int y, float amount);
    void addVelocity(int x, int y, float amountX, float amountY);

    void applyBoundaryConditions(int b, Grid& x) const;
    void linearSolve(int b, Grid& x, const Grid& x0, float a, float c) const;
    void projectVelocity(Grid& velocityX, Grid& velocityY,
                 Grid& p, Grid& div) const;
    void advanceSimulation();

    void advectQuantity(int b, Grid& d, const Grid& d0,
                const Grid& velocityX, const Grid& velocityY, float dt) const;

    void diffuseQuantity(int b, Grid& x, const Grid& x0,
                 float diff, float dt) const;
    void renderDensity(sf::RenderWindow& window) const;
    void decayDensity();
};

