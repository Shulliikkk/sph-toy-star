#include <random>
#include "Processing.h"

Processing::Processing(int N_part, double h, double d, double k, double n, double nu, double lambda, double dt, double maxt) : 
            N_part(N_part), h(h), d(d), k(k), n(n), nu(nu),
            lambda(lambda), dt(dt), maxt(maxt),
            time_steps(maxt / dt), N_grid_cells(d / 2 / h),
            system_radius(std::sqrt(N_grid_cells) * h) {
}

Processing::Processing(int N_part) :
    Processing(N_part, 0.1, 8, 0.1, 1, 1, 4, 0.04, 6) {
}


void Processing::calc_density_for_all(Vec_1d<Particle>& particles, 
                                      Vec_1d<Particle>* grid) {
    for (Particle& particle : particles) {
        particle.change_density(particle.get_m() * particle.W(particle, h));

        int n_grid = N_grid_cells / 2 + particle.get_position()[0] /2/h; 
        int m_grid = N_grid_cells / 2 - particle.get_position()[1] /2/h;

        for (int k = m_grid - 1; k < m_grid + 2; k++) {
            if (k < 0 || k > N_grid_cells - 1) continue;
            for (int l = n_grid - 1; l < n_grid + 2; l++) {
                if (l < 0 || l > N_grid_cells - 1) continue;
				for (Particle grid_particle : grid[N_grid_cells * k + l]) {
					double rho_ij = particle.get_m() * particle.W(grid_particle, h);
					particle.change_density(particle.get_density() + rho_ij);
					grid_particle.change_density(grid_particle.get_density() + rho_ij);
				}
            }
        }
    }
}

void Processing::calc_pressure_for_all(Vec_1d<Particle>& particles) {
    for (Particle& particle : particles) {
        particle.change_pressure(
                k * std::pow(particle.get_density(), 1 + 1 / n)
        );
    }
}

Vec_2d<double> Processing::calc_acceleration_for_all(
        Vec_1d<Particle>& particles,
        Vec_1d<Particle>* grid) {
    calc_density_for_all(particles, grid);
    calc_pressure_for_all(particles);

    Vec_1d<double> curr_density(N_part);
    Vec_2d<double> a(N_part, Vec_1d<double>(2));
    Vec_1d<double> a_p(2);

    for (int i = 0; i < N_part; i++) {
        auto vel = particles[i].get_velocity();
        auto pos = particles[i].get_position();

        a[i][0] += - nu * vel[0] - lambda * pos[0];
        a[i][1] += - nu * vel[1] - lambda * pos[1];

        for (int j = i + 1; j < N_part; j++) {
            auto pres = particles[i].get_pressure();
            auto dens = particles[i].get_density();
            auto j_pres = particles[j].get_pressure();
            auto j_dens = particles[j].get_density();

            auto a_p_coef = - particles[i].get_m() * (pres / (dens * dens)
                            + j_pres / (j_dens * j_dens));
            auto grad = particles[i].gradW(particles[j], h);

            a_p[0] = a_p_coef * grad[0];
            a_p[1] = a_p_coef * grad[1];

            a[i][0] += a_p[0];
            a[i][1] += a_p[1];
            a[j][0] += -a_p[0];
            a[j][1] += -a_p[1];
        }
        curr_density[i] = particles[i].get_density();
    }
    density.push_back(curr_density);
    return a;
}

Vec_2d<Particle> Processing::calc() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> NormRand(0, system_radius);
    
    Vec_1d<Particle> particles;
    for (int i = 0; i < N_part; i++) {
        particles.push_back(Particle(Vec_1d<double>{NormRand(gen),
                                                    NormRand(gen)},
                                     Vec_1d<double>{0, 0},
                                     0.03)
        );
    }

    Vec_2d<Particle> result(time_steps);
    Vec_1d<Particle> grid[N_grid_cells * N_grid_cells];

    for (int i = 0; i < N_part; i++) {
        result[0].push_back(particles[i]);
        int n_grid = N_grid_cells / 2 + particles[i].get_position()[0] / (2 * h);
        int m_grid = N_grid_cells / 2 - particles[i].get_position()[1] / (2 * h);
        grid[N_grid_cells * m_grid + n_grid].push_back(particles[i]);
    }

    Vec_2d<double> acc = calc_acceleration_for_all(particles, grid);
    for (int i = 1; i < time_steps; i++) {
        for (int j = 0; j < N_part; j++) {
            particles[j].change_velocity(acc[j], dt);
            particles[j].change_position(dt);
        }
        acc = calc_acceleration_for_all(particles, grid);

        for (int k = 0; k < N_grid_cells * N_grid_cells; k++) {
            grid[k].clear();
        }

        for (int j = 0; j < N_part; j++) {
            result[i].push_back(particles[j]);
            int n_grid = N_grid_cells / 2 + particles[j].get_position()[0] / (2 * h);
            int m_grid = N_grid_cells / 2 - particles[j].get_position()[1] / (2 * h);

            grid[N_grid_cells * m_grid + n_grid].push_back(particles[j]);
        }
    }

    return result;
}

int Processing::get_N() const { return N_part; }
double Processing::get_system_radius() const { return system_radius; }
Vec_2d<double> Processing::get_density_for_all() { return density; }

