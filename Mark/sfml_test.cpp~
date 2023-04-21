#include <SFML/Graphics.hpp>
#include <iostream> 
#include <string>
#include <typeinfo>
using namespace sf;

template <typename T>
std::ostream& operator << (std::ostream &os, sf::Vector2<T> vec)
{
    return os << "(" << vec.x << ";" << vec.y << ")";
}

class Active_Box {
private:
    sf::Vector2f _position;
    sf::Vector2f _size;
    sf::Color _color;
    std::string _text;
    sf::RenderWindow& _window;
    sf::RectangleShape _rect;
    sf::Font _font;
    sf::Text _text_in_box;
    bool active = false;
    const int size_of_name_text = 40;
    const sf::Vector2f text_offset = sf::Vector2f(50, 25);
    const sf::Color color_of_text = sf::Color::Green;
public:
    Active_Box(sf::RenderWindow& window, 
        sf::Vector2f position, 
        sf::Vector2f size, 
        std::string text_, 
        sf::Color color,
        std::string font_file): 
    _window(window),
    _position(position), 
    _size(size), 
    _text(text_), 
    _color(color)
    {
        sf::RectangleShape rect(_size);
        _rect = rect;
        _rect.setFillColor(_color);
        _rect.setPosition(_position);

        sf::Font font;
        font.loadFromFile(font_file);
        _font = font;

        sf::Text text(_text, _font, size_of_name_text);
        _text_in_box = text;
        _text_in_box.setColor(color_of_text);
        _text_in_box.setPosition(_position + text_offset);
    }

    sf::Vector2f get_position() {
        return _position;
    }

    void show() {
        _window.draw(_rect);
        _window.draw(_text_in_box);
    }
    sf::Vector2f get_size() {
        return _size;
    }
    void set_text(std::string new_text) {
        _text = new_text;
    }
    void set_active() {
        active = true;
    }
    void set_unactive() {
        active = false;
    }
    bool get_active() {
        return active;
    }


};

class Textbox {
private:
    sf::Vector2f _position;
    sf::Vector2f _size;
    std::string _text;
    std::string _name_of_textbox;
    sf::RenderWindow& _window;
    sf::RectangleShape _rect;
    sf::Text _name;
    sf::Font _font;
    sf::Color _color;
    sf::Text _text_in_box;
    const int size_of_name_text = 20;
    const sf::Color color_of_text = sf::Color::Green;
    const sf::Vector2f text_offset = sf::Vector2f(0, -25);
    bool active = false;
public:
    Textbox(sf::RenderWindow& window, 
        sf::Vector2f position, 
        sf::Vector2f size, 
        std::string text_, 
        std::string name_of_textbox, 
        sf::Color color,
        std::string font_file): 
    _window(window),
    _position(position), 
    _size(size), 
    _text(text_), 
    _name_of_textbox(name_of_textbox), 
    _color(color)
    {
        sf::RectangleShape rect(_size);
        _rect = rect;
        _rect.setFillColor(_color);
        _rect.setPosition(_position);

        sf::Font font;
        font.loadFromFile(font_file);
        _font = font;

        sf::Text name(_name_of_textbox, _font, size_of_name_text);
        _name = name;
        _name.setColor(color_of_text);
        _name.setPosition(_position + text_offset);

        sf::Text text(_text, _font, size_of_name_text);
        _text_in_box = text;
        _text_in_box.setColor(color_of_text);
        _text_in_box.setPosition(_position);
    }

    std::string get_text() {
        return _text;
    }
    void add_text(char new_symbol) {
        _text.push_back(new_symbol);
    }

    void delete_text() {
        _text.pop_back();
    }
    
    void update_text() {
        _text_in_box.setString(_text);
    }
    void show() {
        _window.draw(_rect);
        _window.draw(_name);
        _window.draw(_text_in_box);
    }
    sf::Vector2f get_size() {
        return _size;
    }

    void set_text(std::string new_string) {
        _text = new_string;
    }

    sf::Vector2f get_position() {
        return _position;
    }
    void set_active() {
        active = true;
    }
    void set_unactive() {
        active = false;
    }
    bool get_active() {
        return active;
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
                // std::cout << "im here" << std::endl;
                box.set_active();
                box.set_text(" ");
            } 

            if (mouse_position_in_f.x < box.get_position().x || mouse_position_in_f.x > box.get_position().x + box.get_size().x ||
        mouse_position_in_f.y < box.get_position().y || mouse_position_in_f.y > box.get_position().y + box.get_size().y) {
                //std::cout << "im not here" << std::endl;
                box.set_unactive();
            }

        }

}


int main()
{
    sf::Color color = sf::Color::Red;

    sf::Vector2f position(100, 50);
    sf::Vector2f size(400, 30);
    std::string text = "enter the number of particles";
    std::string name = "number of particles:";

    sf::Vector2f position2(700, 50);
    sf::Vector2f size2(400, 30);
    std::string text2 = "enter the time of simulation";
    std::string name2 = "time of simulation:";

    sf::Vector2f position3(250, 300);
    sf::Vector2f size3(400, 30);
    std::string text3 = "enter the radius of particle";
    std::string name3 = "radius:";

    sf::Vector2f position4(700, 300);
    sf::Vector2f size4(400, 80);
    std::string text4 = "click_here_to_start";

    sf::RenderWindow window(sf::VideoMode(1376, 768), "SPH");
    Textbox textbox(window, position, size, text, name, color, "arial.ttf");
    Textbox textbox2(window, position2, size2, text2, name2, color, "arial.ttf");
    Textbox textbox3(window, position3, size3, text3, name3, color, "arial.ttf");
    Active_Box button(window, position4, size4, text4, color, "arial.ttf");

    while (window.isOpen())
    {
 
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
            if (event.type == Event::MouseButtonPressed) {
                check_mouse_click<Textbox>(event, textbox, window);
                check_mouse_click<Textbox>(event, textbox2, window);
                check_mouse_click<Textbox>(event, textbox3, window);
                check_mouse_click<Active_Box>(event, button, window);
            }
            if (textbox.get_active()) {
                if (event.type == sf::Event::TextEntered) {
                    if (event.text.unicode < 128) 
                        textbox.add_text(static_cast<char>(event.text.unicode));
                    textbox.update_text();
                }
            }
            if (textbox2.get_active()) {
                if (event.type == sf::Event::TextEntered) {
                    if (event.text.unicode < 128) 
                        textbox2.add_text(static_cast<char>(event.text.unicode));
                    textbox2.update_text();
                }
            }
            if (textbox3.get_active()) {
                if (event.type == sf::Event::TextEntered) {
                    if (event.text.unicode < 128) 
                        textbox3.add_text(static_cast<char>(event.text.unicode));
                    textbox3.update_text();
                }
            }
            if (button.get_active()) {
                window.close();
            }
                
       window.clear(sf::Color::Blue);
       textbox.show();
       textbox2.show();
       textbox3.show();
       button.show();
       window.display();
       }
    }
    std::cout << "number of particles: " << textbox.get_text() << std::endl;
    std::cout << "time of simulation: " << textbox2.get_text() << std::endl;
    std::cout << "radius of praticle: " << textbox3.get_text() << std::endl;   
    return 0;
}