#include "Processing.h"

void calc_density_for_all(Vec_1d<Particle>& particles,
                          Vec_1d<Particle>* grid,
                          double h, int N_grid_cells) {
    int N = particles.size();

    for (Particle& particle : particles) {
        particle.change_density(particle.get_m() * 
                            particle.W(particle, h)
        );
        int n_grid = N_grid_cells / 2 + particle.get_position()[0] /2/h;
        int m_grid = N_grid_cells / 2 - particle.get_position()[1] /2/h;

        for (int k = m_grid - 1; k < m_grid + 2; k++) {
            if (k < 0 || k > N_grid_cells - 1) continue;
            for (int l = n_grid - 1; l < n_grid + 2; l++) {
                if (l < 0 || l > N_grid_cells - 1) continue;
                //здесь бы функцию
                for (auto grid_particle : grid[N_grid_cells * k + l]) {
                    double rho_ij = grid_particle.get_m() * 
                                    particle.W(grid_particle, h);
                    particle.change_density(particle.get_density() * rho_ij);
                    grid_particle.change_density(grid_particle.get_density() * rho_ij
                    );
                }
            }
        }
    }
}

void calc_pressure_for_all(Vec_1d<Particle>& particles,
                           double k, double n) {
    for (Particle& particle : particles) {
        particle.change_pressure(
                k * std::pow(particle.get_density(), 1 + 1. / n)
        );
    }
}

Vec_2d<double> acceleration(Vec_1d<Particle> particles,
                            Vec_1d<Particle>* grid,
                            double h, double k, double n, double lambda,
                            double nu, double N_grid_cells) {
    int N = particles.size();

    calc_density_for_all(particles, grid, h, N_grid_cells);
    calc_pressure_for_all(particles, k, n);
    Vec_2d<double> a(N, Vec_1d<double>(2));
    Vec_1d<double> a_p(2);

    for (int i = 0; i < N; i++) {
        a[i][0] += - nu * particles[i].get_velocity()[0]
                   - lambda * particles[i].get_position()[0];
        a[i][1] += - nu * particles[i].get_velocity()[1]
                   - lambda * particles[i].get_position()[1];

        for (int j = i + 1; j < N; j++) {
            auto i_pres = particles[i].get_pressure();
            auto i_dens = particles[i].get_density();
            auto j_pres = particles[j].get_pressure();
            auto j_dens = particles[j].get_density();

            auto a_p_coef = - particles[i].get_m() * (i_pres / (i_dens * i_dens) + j_pres / (j_dens * j_dens));
            auto grad = particles[i].gradW(particles[j], h);

            a_p[0] = a_p_coef * grad[0];
            a_p[1] = a_p_coef * grad[1];

            a[i][0] += a_p[0];
            a[i][1] += a_p[1];
            a[j][0] += -a_p[0];
            a[j][1] += -a_p[1];
        }
    }

    return a;
}

Vec_2d<Particle> calc(Vec_1d<Particle>& particles, double h, double d, double k, double n, double lmbda, double nu, double tmax, double dt) {
    int N = particles.size();
    int time_steps = tmax / dt;
    Vec_2d<Particle> result(time_steps);

    const int N_grid_cells = d / (2 * h);
    Vec_1d<Particle> grid[N_grid_cells * N_grid_cells];

    for (Particle& particle : particles) {
        result[0].push_back(particle);
        int n_grid = N_grid_cells / 2 + particle.get_position()[0] / (2 * h);
        int m_grid = N_grid_cells / 2 - particle.get_position()[1] / (2 * h);
        grid[N_grid_cells * m_grid + n_grid].push_back(particle);
    }

    Vec_2d<double> acc = acceleration(particles, grid, h, k, n, lmbda, nu, N_grid_cells);
    Vec_1d<double> rho(N);
    for (int i = 0; i < time_steps; i ++){
        for (int j = 0; j < N; j++){
            particles[j].change_velocity(acc[j], dt);
            particles[j].change_position(particles[j].get_velocity(), dt);
            acc = acceleration(particles, grid, h, k, n, lmbda, nu, N_grid_cells);
            particles[j].change_velocity(acc[j], dt);
        }

        for (int k = 0; k < N_grid_cells * N_grid_cells; k++) {
            grid[k].clear();
        }

        for (Particle& particle : particles) {
            result[i].push_back(particle);
            int n_grid = N_grid_cells / 2 + particle.get_position()[0] / (2 * h);
            int m_grid = N_grid_cells / 2 - particle.get_position()[1] / (2 * h);
            grid[N_grid_cells * m_grid + n_grid].push_back(particle);
        }
    }

    return result;
}

