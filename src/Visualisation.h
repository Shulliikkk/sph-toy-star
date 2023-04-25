#ifndef VISUALISATION_H
#define VISUALISATION_H

#include <SFML/Graphics.hpp>
#include "Particle.h"
#include "Utils.h"

class Visualisation {
private:
    int N_part;
    double h, system_radius;

    Vec_2d<Particle> result;
    Vec_2d<double> density;    
    std::vector<std::pair<double, double>> rho_steps_bounds;

    std::vector<sf::CircleShape> sprites;
    double min_rho, max_rho, rhos_quarter, rhos_half, rhos_step;

    sf::VertexArray quad1;
    sf::VertexArray quad2;
    sf::VertexArray quad3;
    
    sf::Font font;
    sf::Text text_min;
    sf::Text text_max;
    sf::Text text_half;
    sf::Text text_quarter;
    sf::Text text_density;

    Vec_2d<sf::Vertex> lines;

    sf::Color hsv(int hue, double sat, double val);
    void round_string(std::string& string);

public:
    Visualisation(Vec_2d<Particle> result, Vec_2d<double> density,
                  double h, double system_radius);
    void loop();
};

#endif // VISUALISATION_H
