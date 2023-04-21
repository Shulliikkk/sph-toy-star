#include <string>
#include "Textbox.h"

Textbox::Textbox(sf::RenderWindow& window,
    sf::Vector2f position,
    sf::Vector2f size,
    std::string text_,
    std::string name_of_textbox,
    std::string font_file,
    std::string texture_file):
Active_Box(window, position, size, texture_file),
_text(text_),
_name_of_textbox(name_of_textbox)
{
    sf::Font font;
    font.loadFromFile(font_file);
    _font = font;

    sf::Text text(_text, _font, size_of_name_text);
    _text_in_box = text;
    _text_in_box.setColor(color_of_text);
    _text_in_box.setPosition(_position + text_offset);

    sf::Text name(_name_of_textbox, _font, size_of_name_text);
    _name = name;
    _name.setColor(color_of_text);
    _name.setPosition(_position + sf::Vector2f(15, 33));
}

std::string Textbox::get_text() {
    return _text;
}

sf::Vector2f Textbox::get_position() {
    return _position + text_offset;
}

void Textbox::set_text(std::string new_text) {
    _text = new_text;
}

void Textbox::change_text(sf::Event& event) {
    if (active) {
        if (event.type == sf::Event::TextEntered) {
            if (event.text.unicode <= 57 && event.text.unicode >= 48) {
                _text.push_back(static_cast<char>(event.text.unicode));
            }
            else if (event.text.unicode == 8) {
                _text = _text.substr(0, _text.size() - 1);

void Textbox::show() {
    _window.draw(_textbox_sprite);
    _window.draw(_name);
    _window.draw(_text_in_box);
}

void Textbox::text_disappearence() {
    if (active) {
        _text_in_box.setString("");
    }

    else if (active == false && _text.size() == 0) {
        _text_in_box.setString("enter some value");
    }
}
