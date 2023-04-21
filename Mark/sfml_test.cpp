#include <SFML/Graphics.hpp>
#include <iostream> 
#include <string>
#include <cmath>
using namespace sf;


template <typename T>
std::ostream& operator << (std::ostream &os, sf::Vector2<T> vec)
{
    return os << "(" << vec.x << ";" << vec.y << ")";
}


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

    sf::Vector2f get_position() {
        return _position;
    }

    void show() {
        _window.draw(_textbox_sprite);
    }

    sf::Vector2f get_size() {
        return _size;
    }

    void set_active() {
        active = true;
    }
    void set_text(std::string new_text) {
        _text = new_text;
    }
    void set_unactive() {
        active = false;
    }

    bool get_active() {
        return active;
    }
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

    std::string get_text() {
        return _text;
    }

    sf::Vector2f get_position() {
        return _position + text_offset;
    }

    void set_text(std::string new_text) {
        _text = new_text;
    }

    void change_text(sf::Event& event) {
        if (active) {
            if (event.type == sf::Event::TextEntered) {
                if (event.text.unicode <= 57 && event.text.unicode >= 48) {
                    _text.push_back(static_cast<char>(event.text.unicode));
                }
                else if (event.text.unicode == 8) {
                    _text = _text.substr(0, _text.size() - 1);
                }
                _text_in_box.setString(_text);
            }  
        }
    }
    
    void show() {
        _window.draw(_textbox_sprite);
        _window.draw(_name);
        _window.draw(_text_in_box);
    }

    void text_disappearence() {
        if (active) {
            _text_in_box.setString("");
        }

        else if (active == false && _text.size() == 0) {
            _text_in_box.setString("enter some value");
        }
    }
};

template <typename T>
bool check_mouse_click(sf::Event& event, T& box, sf::RenderWindow& window) {
    sf::Vector2f mouse_position_in_f;
    sf::Vector2i mouse_position = sf::Mouse::getPosition(window);
    mouse_position_in_f.x = mouse_position.x;
    mouse_position_in_f.y = mouse_position.y; 
    if (event.key.code == Mouse::Left) {
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

int string2int(std::string string_) {
    int size = string_.size();
    int result_number = 0;
    for (auto i = 0; i < size; i++) {
        result_number += std::pow(10, size - i - 1) * (static_cast<int>(string_[i]) - 48);
    }
    return result_number;
}

sf::Vector3i menu()
{   
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
            if (event.type == Event::MouseButtonPressed) {
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
    
    sf::Vector3i result(string2int(textbox.get_text()), string2int(textbox2.get_text()), string2int(textbox3.get_text()));
    return result;
}


int main() {
    sf::Vector3i menu_data;

    menu_data = menu();
    std::cout << menu_data.x << ' ' << menu_data.y << ' ' << menu_data.z << std::endl;
}