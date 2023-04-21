#ifndef PARTICLE_H
#define PARTICLE_H

#include "Utils.h"

class Particle {
private:
    Vec_1d<double> velocity = Vec_1d<double>(3);
    Vec_1d<double> position = Vec_1d<double>(3);
    double m;
    double rho;
    double p;

public:
    Particle(Vec_1d<double> position, Vec_1d<double> velocity, double m);

    double distance_to(Particle& other) const;
    double W(Particle& other, double h);
    Vec_1d<double> gradW(Particle& other, double h);

    Vec_1d<double> get_position() const;
    Vec_1d<double> get_velocity() const;
    double get_m() const;
    double get_density() const;
    double get_pressure() const;

    void change_velocity(Vec_1d<double>& acceleration, double dt);
    void change_position(Vec_1d<double> velocity, double dt);
    void change_density(double new_density);
    void change_pressure(double new_pressure);
};

#endif // PARTICLE_H
