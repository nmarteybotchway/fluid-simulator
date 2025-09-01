#include <SFML/Graphics.hpp>
#include "Fluid.h"
#include <chrono>
#include <iostream>

int main() {
    int n = 256;
    Fluid fluid(n, 4, 0.1, 0, 0, 3);
    sf::RenderWindow window(sf::VideoMode(n * fluid.getScale(), n * fluid.getScale()), "Fluid Simulation");

    sf::Vector2i prevMousePos = sf::Mouse::getPosition(window);
    bool firstFrame = true;

    double totalTime = 0.0;
    int stepCount = 0;

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        sf::Vector2i mousePos = sf::Mouse::getPosition(window);

        float amtX = 0.0f;
        float amtY = 0.0f;
        if (!firstFrame) {
            amtX = static_cast<float>(mousePos.x - prevMousePos.x);
            amtY = static_cast<float>(mousePos.y - prevMousePos.y);
        } else {
            firstFrame = false;
        }
        prevMousePos = mousePos;

        if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
            int xCell = mousePos.x / fluid.getScale();
            int yCell = mousePos.y / fluid.getScale();

            if (xCell >= 0 && xCell < n && yCell >= 0 && yCell < n) {
                fluid.addDensity(xCell, yCell, 100);
                fluid.addVelocity(xCell, yCell, amtX / fluid.getScale(), amtY / fluid.getScale());
            }
        }

        window.clear();
        fluid.renderD(window);
        fluid.fadeD();
        window.display();
    }
}
