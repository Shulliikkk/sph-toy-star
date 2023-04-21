#ifndef PROCESSING_H
#define PROCESSING_H

#include <iostream>
#include "Utils.h"
#include "Particle.h"

class Processing {
private:
    int N_part;
    double h, d, k, n, nu, lambda, dt, maxt;
    const int time_steps, N_grid_cells;
    Vec_2d<double> density;

public:
    Processing(int N_part); // default parameters
    Processing(int N_part, double h, double d, double k, double n, double nu, double lambda, double dt, double maxt);

    void calc_density_for_all(Vec_1d<Particle>& particles, Vec_1d<Particle>* grid);
    void calc_pressure_for_all(Vec_1d<Particle>& particles);
    Vec_2d<double> calc_acceleration_for_all(Vec_1d<Particle>& particles, Vec_1d<Particle>* grid);
    Vec_2d<Particle> calc();

    int get_N() const;
    double get_dt() const;
    double get_maxt() const;
    Vec_2d<double> get_density_for_all();
};

#endif // PROCESSING_H
