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

int main() {
    Processing model(100, 0.1, 20, 0.1, 1, 1, 4, 0.04, 12);

    Vec_2d<Particle> result;
    {
        Timer<std::chrono::seconds> t;
        result = model.calc();
    }

    Visualisation visualisation;
    Vec_2d<double> density = model.get_density_for_all();
    visualisation.loop(result, density, 0.1);
    
    return 0;
}
