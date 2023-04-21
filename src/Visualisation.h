#ifndef VISUALISATION_H
#define VISUALISATION_H

#include <SFML/Graphics.hpp>
#include "Particle.h"
#include "Utils.h"

class Visualisation {
private:
    sf::Color hsv(int hue, double sat, double val);

public:
    Visualisation();
    void loop(Vec_2d<Particle> result, Vec_2d<double> density, double h);
};

#endif // VISUALISATION_H
