#include <SFML/Graphics.hpp>
#include "Fluid.h"

int main() {
    int n = 128;
    int scale = 10;
    Fluid fluid(n, 4, 0.1, 0, 0);
    sf::RenderWindow window(sf::VideoMode(n * scale, n * scale), "Fluid Simulation");

    sf::Vector2i prevMousePos = sf::Mouse::getPosition(window);
    bool firstFrame = true;

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        // Get current mouse position each frame
        sf::Vector2i mousePos = sf::Mouse::getPosition(window);

        // Compute velocity delta
        float amtX = 0.0f;
        float amtY = 0.0f;
        if (!firstFrame) {
            amtX = static_cast<float>(mousePos.x - prevMousePos.x);
            amtY = static_cast<float>(mousePos.y - prevMousePos.y);
        } else {
            firstFrame = false;
        }

        prevMousePos = mousePos; // update for next frame

        // If left mouse pressed, apply to fluid
        if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
            int xCell = mousePos.x / scale;
            int yCell = mousePos.y / scale;
            fluid.addDensity(xCell, yCell, 100);
            fluid.addVelocity(xCell, yCell, amtX, amtY);
        }

        fluid.step();
        window.clear();
        fluid.renderD(window);
        fluid.fadeD();
        window.display();
    }
}
