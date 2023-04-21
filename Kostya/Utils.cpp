#include "Utils.h"

template<typename T>
std::ostream& operator << (std::ostream &os, sf::Vector2<T> vec)
{
    return os << "(" << vec.x << ";" << vec.y << ")";
}
