#include <SFML/Graphics.hpp>
#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <random>
#include <unistd.h>
#include <chrono>

#include "Utils.h"
#include "Timer.h"
#include "Particle.h"
#include "Processing.h"

int main() {
    Processing model(50);
    //Processing::Processing(int N_part, double h, double d, double k, double n, double nu, double lambda, double dt, double maxt) : 

    Vec_2d<Particle> result;
    Vec_1d<Particle> particles;
    Vec_2d<double> a(model.get_N(), Vec_1d<double>(2));

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> NormRand(0, 0.1);

    {
        Timer<std::chrono::seconds> t;

        for (int i = 0; i < model.get_N(); i++) {
            particles.push_back(Particle(Vec_1d<double>{NormRand(gen),
                                                        NormRand(gen), 0},
                                         Vec_1d<double>{NormRand(gen),
                                                        NormRand(gen), 0},
                                         0.03)
            );
        }

        result = model.calc(particles);
    }

    std::vector<sf::CircleShape> sprites(model.get_N());
    std::size_t step = 0;
    sf::RenderWindow window(sf::VideoMode (400, 400), " SFML works!");
    double want_fps = 5;
    sf::Clock loop_timer;

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }
        window.clear();
        for (sf::CircleShape &sprite : sprites) {
            sprite.setRadius(1);
            for(std::size_t j = 0; j < result[step].size(); j++) {
                sprite.setPosition(result[step][j].get_position()[0] * 100 + 200,
                                   result[step][j].get_position()[1] * 60 + 200);
                window.draw(sprite);
            }
        }

        window.display();

        if (step + 1 < result.size()) step += 1;
        else step = 0;

        sf::Int32 frame_duration = loop_timer.getElapsedTime().asMilliseconds();
        sf::Int32 time_to_sleep = int(1000.f/want_fps) - frame_duration;
        if (time_to_sleep > 0) {
            sf::sleep(sf::milliseconds(time_to_sleep));
        }
        loop_timer.restart();
    }

    return 0;
}
