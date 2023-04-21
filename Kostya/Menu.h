#ifndef MENU_H
#define MENU_H
#include <SFML/Graphics.hpp>
#include <string>
#include "Textbox.h"

class Menu {
private:
  template <typename T>
  bool check_mouse_click(sf::Event& event, T& box, sf::RenderWindow& window);
  int string2int(std::string string_);

public:
  Menu() {};
  sf::Vector3i menu();
};

#endif
