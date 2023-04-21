#ifndef ACTIVE_BOX_H
#define ACTIVE_BOX_H

#include <SFML/Graphics.hpp>
#include "Utils.h"

class Active_Box {
private:
    std::string _text;
protected:
    std::string _texture_file;
    sf::Vector2f _position;
    sf::Vector2f _size;
    sf::RenderWindow& _window;
    sf::Texture _texture;
    sf::Sprite _textbox_sprite;
    bool active = false;
    const sf::Color color_of_text = sf::Color::Black;

public:
    Active_Box(sf::RenderWindow& window,
        sf::Vector2f position,
        sf::Vector2f size,
        std::string texture_file);

    sf::Vector2f get_position();
    sf::Vector2f get_size();
    bool get_active();

    void set_active();
    void set_text(std::string new_text);
    void set_unactive();

    void show();
};

#endif
