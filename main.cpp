#include <iostream>
#include <SFML/Graphics.hpp>
#include "Fluid.h"

int main() {
    int n{64};
    int scale = 10;
    Fluid fluid = Fluid(n, 4, 00.1, 0, 0);
    sf::RenderWindow window(sf::VideoMode(n*scale, n*scale), "Fluid Simulation");
    bool isMousePressed{false};
    while (window.isOpen()) {
        sf::Event event{};
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }

            if (event.type == sf::Event::MouseButtonPressed) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    isMousePressed = true;
                    fluid.addDensity(event.mouseButton.x/scale, event.mouseButton.y/scale, 100);
                }
            }

            if (event.type == sf::Event::MouseButtonReleased) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    isMousePressed = false;
                }
            }
        }

        if (isMousePressed) {
            sf::Vector2i mousePos{sf::Mouse::getPosition(window)};
            fluid.addDensity(mousePos.x/scale, mousePos.y/scale, 100);
        }

        fluid.step();

        window.clear();
        fluid.renderD(window);
        window.display();
    }
}