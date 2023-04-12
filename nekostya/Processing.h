#ifndef PROCESSING_H
#define PROCESSING_H

#include "Utils.h"
#include "Particle.h"

void calc_density_for_all(Vec_1d<Particle>& particles,
                          Vec_1d<Particle>* grid,
                          double h, int N_grid_cells);

void calc_pressure_for_all(Vec_1d<Particle>& particles,
                           double k, double n);

Vec_2d<double> acceleration(Vec_1d<Particle> particles,
                            Vec_1d<Particle>* grid,
                            double h, double k, double n, double lambda,
                            double nu, double N_grid_cells);

Vec_2d<Particle> calc(Vec_1d<Particle>& particles, double h, double d, double k, double n, double lmbda, double nu, double tmax, double dt);

#endif // PROCESSING_H
