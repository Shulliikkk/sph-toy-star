#ifndef TEXTBOX_H
#define TEXTBOX_H

#include <SFML/Graphics.hpp>
#include "Active_Box.h"

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
        std::string texture_file);

    std::string get_text();

    sf::Vector2f get_position();

    void set_text(std::string new_text);

    void change_text(sf::Event& event);

    void show();

    void text_disappearence();
};

#endif
