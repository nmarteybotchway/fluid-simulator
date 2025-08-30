#include <iostream>
#include <SFML/Graphics.hpp>
#include "Fluid.h"

int main() {
    int n = 64;
    int scale = 10;
    Fluid fluid(n, 4, 0.1, 0, 0);
    sf::RenderWindow window(sf::VideoMode(n * scale, n * scale), "Fluid Simulation");

    sf::Vector2i prevMousePos = sf::Mouse::getPosition(window);
    bool firstFrame = true;
    bool isMousePressed = false;

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();

            if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left)
                isMousePressed = true;

            if (event.type == sf::Event::MouseButtonReleased && event.mouseButton.button == sf::Mouse::Left)
                isMousePressed = false;
        }

        // Get current mouse position every frame
        sf::Vector2i mousePos = sf::Mouse::getPosition(window);

        // Compute delta (velocity)
        float amtX = 0.0f;
        float amtY = 0.0f;
        if (!firstFrame) {
            amtX = static_cast<float>(mousePos.x - prevMousePos.x);
            amtY = static_cast<float>(mousePos.y - prevMousePos.y);
        } else {
            firstFrame = false;
        }
        prevMousePos = mousePos;

        // Apply to fluid if left mouse is pressed
        if (isMousePressed) {
            fluid.addDensity(mousePos.x / scale, mousePos.y / scale, 100);
            fluid.addVelocity(mousePos.x / scale, mousePos.y / scale, amtX, amtY);
            std::cout << amtX <<amtY <<std::endl;
        }

        // Step the fluid simulation
        fluid.step();

        // Render
        window.clear();
        fluid.renderD(window);
        window.display();
    }
}
