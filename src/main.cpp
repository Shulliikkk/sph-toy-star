#include <SFML/Graphics.hpp>
#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <unistd.h>
#include <chrono>

#include "Utils.h"
#include "Timer.h"
#include "Particle.h"
#include "Processing.h"
#include "Visualisation.h"
#include "Menu.h"

int main() {

    sf::Vector3i menu_data = call_menu();
    int N_part = menu_data.x;
    double maxt = (double) menu_data.y;
    double h = (double) menu_data.z / 100.;

    Processing model(N_part, h, 20, 0.1, 1, 1, 4, 0.04, maxt);
    //Processing::Processing(int N_part, double h, double d, double k, double n, double nu, double lambda, double dt, double maxt)

    Vec_2d<Particle> result;
    {
        //Timer<std::chrono::seconds> t;
        result = model.calc();
    }

    Visualisation visualisation;
    Vec_2d<double> density = model.get_density_for_all();
    visualisation.loop(result, density, 0.1);
    
    return 0;
}
