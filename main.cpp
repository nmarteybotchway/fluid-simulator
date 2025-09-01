#include <SFML/Graphics.hpp>
#include "Fluid.h"

// Helper function to apply mouse input to the fluid
void applyMouseInteraction(Fluid& fluid, const sf::Vector2i& mousePos, const sf::Vector2i& delta) {
    if (sf::Vector2i cellPos(mousePos.x / fluid.getScale(), mousePos.y / fluid.getScale()); cellPos.x >= 0 && cellPos.x < fluid.getGridSize() &&
                                                                                            cellPos.y >= 0 && cellPos.y < fluid.getGridSize()) {
        fluid.addDensity(cellPos.x, cellPos.y, 100.0f);
        fluid.addVelocity(cellPos.x, cellPos.y,
                          static_cast<float>(delta.x) / fluid.getScale(),
                          static_cast<float>(delta.y) / fluid.getScale());
        }
}

int main() {
    constexpr int gridSize = 256;
    constexpr int fluidScale = 3;
    Fluid fluid(gridSize, 4, 0.1f, 0.0f, 0.0f, fluidScale);

    sf::RenderWindow window(sf::VideoMode(fluid.getWindowWidth(), fluid.getWindowHeight()),
                            "Fluid Simulation",
                            sf::Style::Titlebar | sf::Style::Close);

    sf::Vector2i prevMousePos = sf::Mouse::getPosition(window);

    while (window.isOpen()) {
        sf::Event event{};
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        // Mouse movement delta
        sf::Vector2i mousePos = sf::Mouse::getPosition(window);
        sf::Vector2i delta = mousePos - prevMousePos;
        prevMousePos = mousePos;

        // Apply interaction if left button is pressed
        if (sf::Mouse::isButtonPressed(sf::Mouse::Left))
            applyMouseInteraction(fluid, mousePos, delta);

        // Advance simulation
        fluid.advanceSimulation();

        // Render
        window.clear();
        fluid.renderDensity(window);
        fluid.decayDensity();
        window.display();
    }
}
