#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <vector>

template<typename T>
using Vec_1d = std::vector<T>;
template<typename T>
using Vec_2d = std::vector<std::vector<T>>;
template<typename T>
using Vec_3d = std::vector<std::vector<std::vector<T>>>;

const double PI = std::acos(-1.0);

/*
class Vec_1d {
private:
    double x, y, z;

public:
    Vector3();
    Vector3(double x_coord, double y_coord, double z_coord);
    double length() const;

    Vector3 operator*(double scalar);
    Vector3 operator/(double scalar);
    Vector3 operator+(Vector3 another);
    Vector3 operator-(Vector3 another);
    double operator*(Vector3 another);
    
    void normalise();
    Vector3 normalised();
    
    void print() const;
};
*/

#endif // UTILS_H
