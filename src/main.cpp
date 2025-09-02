#include <SFML/Graphics.hpp>
#include "../include/Fluid.h"

/**
 * @brief Apply mouse interaction to the fluid simulation.
 *
 * Converts the mouse screen coordinates to fluid grid coordinates, then
 * adds density and velocity based on the mouse position and movement.
 *
 * @param fluid Reference to the Fluid simulation object.
 * @param mousePos Current mouse position in window coordinates.
 * @param delta Mouse movement delta since the last frame.
 */
void applyMouseInteraction(Fluid& fluid, const sf::Vector2i& mousePos, const sf::Vector2i& delta) {
    // Convert mouse position to fluid grid coordinates
    if (sf::Vector2i cellPos(mousePos.x / fluid.getScale(), mousePos.y / fluid.getScale());
        cellPos.x >= 0 && cellPos.x < fluid.getGridSize() && cellPos.y >= 0 && cellPos.y < fluid.getGridSize()) {

        // Add density at this grid cell
        fluid.addDensity(cellPos.x, cellPos.y, 100.0f);

        // Add velocity based on mouse movement
        fluid.addVelocity(cellPos.x, cellPos.y,
                          static_cast<float>(delta.x) / fluid.getScale(),
                          static_cast<float>(delta.y) / fluid.getScale());
    }
}

/**
 * @brief Entry point for the fluid simulation program.
 *
 * Sets up an SFML window, initializes the fluid simulation, and runs
 * the main loop handling input, simulation updates, and rendering.
 *
 * @return int Exit status code (0 for success)
 */
int main() {
    constexpr int gridSize = 256;    ///< Number of cells per dimension
    constexpr int fluidScale = 3;    ///< Display scaling factor
    Fluid fluid(gridSize, 4, 0.1f, 0.0f, 0.0f, fluidScale); ///< Initialize fluid simulation

    // Create SFML window
    sf::RenderWindow window(sf::VideoMode(fluid.getWindowWidth(), fluid.getWindowHeight()),
                            "Fluid Simulation",
                            sf::Style::Titlebar | sf::Style::Close);

    // Previous mouse position to compute movement delta
    sf::Vector2i prevMousePos = sf::Mouse::getPosition(window);

    while (window.isOpen()) {
        sf::Event event{};
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }


        // Compute mouse movement delta
        sf::Vector2i mousePos = sf::Mouse::getPosition(window);
        sf::Vector2i delta = mousePos - prevMousePos;
        prevMousePos = mousePos;

        // Apply interaction if left mouse button is pressed
        if (sf::Mouse::isButtonPressed(sf::Mouse::Left))
            applyMouseInteraction(fluid, mousePos, delta);

        // Advance fluid simulation one timestep
        fluid.advanceSimulation();

        // Render density field
        window.clear();
        fluid.renderDensity(window);

        // Optional: decay density over time
        // fluid.decayDensity();

        window.display();
    }

    return 0;
}
