#ifndef UTILS_H
#define UTILS_H

#include <SFML/Graphics.hpp>
#include <cmath>
#include <vector>

template<typename T>
using Vec_1d = std::vector<T>;
template<typename T>
using Vec_2d = std::vector<std::vector<T>>;
template<typename T>
using Vec_3d = std::vector<std::vector<std::vector<T>>>;
template<typename T>
std::ostream& operator << (std::ostream &os, sf::Vector2<T> vec);


const double PI = std::acos(-1.0);

#endif
