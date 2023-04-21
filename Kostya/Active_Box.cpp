#include "Active_Box.h"

Active_Box::Active_Box(sf::RenderWindow& window,
    sf::Vector2f position,
    sf::Vector2f size,
    std::string texture_file):
_window(window),
_position(position),
_size(size),
_texture_file(texture_file)
{
    sf::Texture textbox_texture;
    textbox_texture.loadFromFile(texture_file);
    _texture = textbox_texture;

    _textbox_sprite.setTexture(_texture);
    _textbox_sprite.setPosition(_position);
}

sf::Vector2f Active_Box::get_position() {
    return _position;
}

void Active_Box::show() {
    _window.draw(_textbox_sprite);
}

sf::Vector2f Active_Box::get_size() {
    return _size;
}

void Active_Box::set_active() {
    active = true;
}
void Active_Box::set_text(std::string new_text) {
    _text = new_text;
}
void Active_Box::set_unactive() {
    active = false;
}

bool Active_Box::get_active() {
    return active;
}
