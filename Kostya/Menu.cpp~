#include <string>
#include <cmath>
#include <SFML/Graphics.hpp>
#include "Textbox.h"
#include "Active_Box.h"
#include "Menu.h"

template <typename T>
bool Menu::check_mouse_click(sf::Event& event, T& box, sf::RenderWindow& window) {
    sf::Vector2f mouse_position_in_f;
    sf::Vector2i mouse_position = sf::Mouse::getPosition(window);
    mouse_position_in_f.x = mouse_position.x;
    mouse_position_in_f.y = mouse_position.y;
    if (event.key.code == sf::Mouse::Left) {
        if (mouse_position_in_f.x >= box.get_position().x && mouse_position_in_f.x <= box.get_position().x + box.get_size().x &&
            mouse_position_in_f.y >= box.get_position().y && mouse_position_in_f.y <= box.get_position().y + box.get_size().y) {
                box.set_text("");
                box.set_active();
        }

        if (mouse_position_in_f.x < box.get_position().x || mouse_position_in_f.x > box.get_position().x + box.get_size().x ||
            mouse_position_in_f.y < box.get_position().y || mouse_position_in_f.y > box.get_position().y + box.get_size().y) {
                box.set_unactive();
        }
    }
}

int Menu::string2int(std::string string_) {
    int size = string_.size();
    int result_number = 0;
    for (auto i = 0; i < size; i++) {
        result_number += std::pow(10, size - i - 1) * (static_cast<int>(string_[i]) - 48);
    }
    return result_number;
}

//Menu::Menu() {};

sf::Vector3i Menu::menu() {
    sf::RenderWindow window(sf::VideoMode(1376, 768), "SPH");

    sf::Color color = sf::Color::Red;

    sf::Vector2f position(100, 50);
    sf::Vector2f size(320, 30);
    std::string text = "enter the number of particles";
    std::string name = "number of particles:";

    std::string texture_file = "textbox_texture3.png";
    std::string texture_file2 = "button_texture.png";

    sf::Vector2f position2(500, 50);
    sf::Vector2f size2(320, 30);
    std::string text2 = "enter the time of simulation";
    std::string name2 = "time of simulation:";

    sf::Vector2f position3(900, 50);
    sf::Vector2f size3(320, 30);
    std::string text3 = "enter the radius of particle";
    std::string name3 = "radius:";

    sf::Vector2f position4(500, 600);
    sf::Vector2f size4(400, 150);
    std::string text4 = "click here to start";

    Textbox textbox(window, position, size, text, name, "arial.ttf", texture_file);
    Textbox textbox2(window, position2, size2, text2, name2, "arial.ttf", texture_file);
    Textbox textbox3(window, position3, size3, text3, name3, "arial.ttf", texture_file);
    Active_Box button(window, position4, size4, texture_file2);

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
            if (event.type == sf::Event::MouseButtonPressed) {
                check_mouse_click<Textbox>(event, textbox, window);
                check_mouse_click<Textbox>(event, textbox2, window);
                check_mouse_click<Textbox>(event, textbox3, window);
                check_mouse_click<Active_Box>(event, button, window);

                if (button.get_active()) {
                    window.close();
                }

                textbox.text_disappearence();
                textbox2.text_disappearence();
                textbox3.text_disappearence();
            }
            textbox.change_text(event);
            textbox2.change_text(event);
            textbox3.change_text(event);


           window.clear(sf::Color::Blue);
           textbox.show();
           textbox2.show();
           textbox3.show();
           button.show();
           window.display();
       }
    }
}

sf::Vector3i Menu::result(string2int(textbox.get_text()), string2int(textbox2.get_text()), string2int(textbox3.get_text()));
    return result;
}
