#ifndef MENU_H
#define MENU_H

#include <SFML/Graphics.hpp>
#include <iostream>
#include <string>
#include <cmath>

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
    bool active;
    const sf::Color color_of_text = sf::Color::Black;

public:
    Active_Box(sf::RenderWindow& window,
               sf::Vector2f position,
               sf::Vector2f size,
               std::string texture_file
    );

    sf::Vector2f get_position();
    void show();
    sf::Vector2f get_size();
    void set_active();
    void set_inactive();
    bool get_active();
    void set_text(std::string new_text);
};

class Textbox : public Active_Box {
private:
    const int size_of_name_text = 20;
    std::string _name_of_textbox;
    sf::Text _name;
    sf::Font _font;
    sf::Text _text_in_box;
    std::string _text;
    const sf::Vector2f text_offset = sf::Vector2f(15, 63);
public:
    Textbox(sf::RenderWindow& window,
            sf::Vector2f position, 
            sf::Vector2f size, 
            std::string text_, 
            std::string name_of_textbox, 
            std::string font_file,
            std::string texture_file
    );

    std::string get_text();
    void set_text(std::string new_text);
    void change_text(sf::Event& event);

    sf::Vector2f get_position();
    void show();
    void text_disappearence();
};

template <typename T>
void check_mouse_click(sf::Event& event,
                       T& box,
                       sf::RenderWindow& window);

int string2int(std::string string_);

sf::Vector3i call_menu();

#endif // MENU_H
