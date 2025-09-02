#pragma once
#include <cmath>
#include <SFML/Graphics.hpp>
#include "Grid.h"

/**
 * @brief 2D fluid simulation using stable fluid solver (Jos Stam style).
 *
 * Stores density and velocity fields on a uniform grid. Supports
 * diffusion, advection, and velocity projection. Optional rendering
 * using SFML textures.
 */
class Fluid {
    int gridSize;              ///< Number of cells per dimension
    int solverIterations;      ///< Number of linear solver iterations per timestep
    float timeStep;            ///< Simulation timestep
    float diffusionRate;       ///< Diffusion rate for density
    float viscosity;           ///< Viscosity for velocity
    int scale;                 ///< Rendering scale factor
    Grid prevDensity;          ///< Temporary density buffer
    Grid density;              ///< Density field
    Grid velocityX;            ///< X component of velocity
    Grid velocityY;            ///< Y component of velocity
    Grid velocityXPrev;        ///< Temporary X velocity buffer
    Grid velocityYPrev;        ///< Temporary Y velocity buffer

    // Reused buffers to avoid per-frame allocation
    mutable sf::Texture densityTex;       ///< Texture for rendering density
    mutable bool densityTexInit = false;  ///< Flag if texture initialized
    mutable std::vector<sf::Uint8> pixels;///< Pixel buffer for texture updates

public:
    /**
     * @brief Construct a new Fluid simulation instance.
     *
     * @param gridSize Number of cells per dimension
     * @param solverIterations Number of solver iterations per timestep
     * @param timeStep Simulation timestep
     * @param diffusionRate Density diffusion rate
     * @param viscosity Fluid viscosity
     * @param scale Rendering scale factor
     */
    Fluid(int gridSize, int solverIterations, float timeStep, float diffusionRate, float viscosity, int scale);

    /// @return The display scale factor
    [[nodiscard]] int getScale() const { return scale; }

    /// @return The number of cells per grid dimension
    [[nodiscard]] int getGridSize() const { return gridSize; }

    /// @return Width of the window in pixels
    [[nodiscard]] int getWindowWidth() const { return gridSize * scale; }

    /// @return Height of the window in pixels
    [[nodiscard]] int getWindowHeight() const { return gridSize * scale; }

    /**
     * @brief Add density at a specific grid cell.
     *
     * @param x X-coordinate of the cell
     * @param y Y-coordinate of the cell
     * @param amount Amount of density to add
     */
    void addDensity(int x, int y, float amount);

    /**
     * @brief Add velocity to a specific grid cell.
     *
     * @param x X-coordinate of the cell
     * @param y Y-coordinate of the cell
     * @param amountX Amount to add to X velocity
     * @param amountY Amount to add to Y velocity
     */
    void addVelocity(int x, int y, float amountX, float amountY);

    /**
     * @brief Apply boundary conditions to a quantity.
     *
     * @param b Type of quantity (0=density, 1=X velocity, 2=Y velocity)
     * @param x Grid to apply boundary conditions to
     */
    void applyBoundaryConditions(int b, Grid& x) const;

    /**
     * @brief Solve a linear system for diffusion or projection.
     *
     * Uses Gauss-Seidel iteration with boundary conditions.
     * @param b Boundary type
     * @param x Grid to solve for
     * @param x0 Right-hand side grid
     * @param a Coefficient for neighbors
     * @param c Diagonal coefficient
     */
    void linearSolve(int b, Grid& x, const Grid& x0, float a, float c) const;

    /**
     * @brief Make the velocity field divergence-free.
     *
     * Projects the velocity to be mass-conserving using pressure projection.
     * @param velocityX X component of velocity
     * @param velocityY Y component of velocity
     * @param p Pressure buffer
     * @param div Divergence buffer
     */
    void projectVelocity(Grid& velocityX, Grid& velocityY, Grid& p, Grid& div) const;

    /**
     * @brief Advance the simulation by one timestep.
     *
     * Performs diffusion, advection, and velocity projection for density and velocity fields.
     */
    void advanceSimulation();

    /**
     * @brief Advect a quantity through the velocity field.
     *
     * @param b Boundary type
     * @param d Destination grid
     * @param d0 Source grid
     * @param velocityX X velocity field
     * @param velocityY Y velocity field
     * @param dt Timestep
     */
    void advectQuantity(int b, Grid& d, const Grid& d0, const Grid& velocityX, const Grid& velocityY, float dt) const;

    /**
     * @brief Diffuse a quantity over time.
     *
     * @param b Boundary type
     * @param x Destination grid
     * @param x0 Source grid
     * @param diff Diffusion coefficient
     * @param dt Timestep
     */

    void diffuseQuantity(int b, Grid& x, const Grid& x0, float diff, float dt) const;

    /**
     * @brief Render the density field to an SFML window.
     *
     * Uses an internal texture and pixel buffer for efficient rendering.
     * @param window SFML render window
     */
    void renderDensity(sf::RenderWindow& window) const;

    /**
     * @brief Decay density values over time.
     *
     * Reduces density slightly per timestep. Typically used for visual effect.
     */
    void decayDensity();
};
