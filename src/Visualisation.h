#ifndef VISUALISATION_H
#define VISUALISATION_H

#include <SFML/Graphics.hpp>
#include "Particle.h"
#include "Utils.h"

class Visualisation {
private:
    int N_part;
    double h, system_radius;
    sf::Color hsv(int hue, double sat, double val);

public:
    Visualisation(double h, double system_radius);
    void loop(Vec_2d<Particle> result, Vec_2d<double> density);
};

#endif // VISUALISATION_H
