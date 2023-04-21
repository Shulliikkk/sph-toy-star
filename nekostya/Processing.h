#ifndef PROCESSING_H
#define PROCESSING_H

#include "Utils.h"
#include "Particle.h"

class Processing {
private:
    const int N_part;
    double h, d, k, n, nu, lambda, dt, maxt;
    const int time_steps, N_grid_cells;

public:
    Processing(int N_part);
    Processing(int N_part, double h, double d, double k, double n, double nu, double lambda, double dt, double maxt);

    void calc_density_for_all(Vec_1d<Particle>& particles,
                              Vec_1d<Particle>* grid);

    void calc_pressure_for_all(Vec_1d<Particle>& particles);

    Vec_2d<double> acceleration(Vec_1d<Particle> particles,
                                Vec_1d<Particle>* grid);

    Vec_2d<Particle> calc(Vec_1d<Particle>& particles);

    int get_N() const;
    double get_dt() const;
    double get_maxt() const;
};

#endif // PROCESSING_H
