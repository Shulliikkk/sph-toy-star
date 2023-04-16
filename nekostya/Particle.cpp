#include "Particle.h"

Particle::Particle(Vec_1d<double> position, Vec_1d<double> velocity, double m) : m(m), velocity(velocity), position(position) {}

double Particle::distance_to(Particle& other) const {
    return std::sqrt(std::pow(position[0] - other.position[0], 2) +
                     std::pow(position[1] - other.position[1], 2) +
                     std::pow(position[2] - other.position[2], 2)
           );
}

double Particle::W(Particle& other, double h) {
    double abs_r = distance_to(other);
    return 1 / std::pow(h, 3) / std::pow(PI, 1.5) *
           std::exp(-abs_r / h * abs_r / h);
}

Vec_1d<double> Particle::gradW(Particle& other, double h) {
    double abs_r = distance_to(other);
    double scalar_deriv = -2 / std::pow(h, 5) / std::pow(PI, 1.5) *
                          std::exp(-abs_r / h * abs_r / h);
    
    Vec_1d<double> grad = {
        scalar_deriv * (position[0] - other.position[0]),
        scalar_deriv * (position[1] - other.position[1]),
        scalar_deriv * (position[2] - other.position[2])
    };
    return grad;
}

Vec_1d<double> Particle::get_position() const { return position; }
Vec_1d<double> Particle::get_velocity() const { return velocity; }
double Particle::get_m() const { return m; }
double Particle::get_density() const { return rho; }
double Particle::get_pressure() const { return p; }

void Particle::change_velocity(Vec_1d<double>& acceleration, double dt) {
    for (int k = 0; k < 3; k++) velocity[k] += acceleration[k] * dt / 2;
}

void Particle::change_position(Vec_1d<double> velocity, double dt) {
    for (int k = 0; k < 3; k++) position[k] += velocity[k] * dt;
}

void Particle::change_density(double new_density) {
    rho = new_density;
}

void Particle::change_pressure(double new_pressure) {
    p = new_pressure;
}