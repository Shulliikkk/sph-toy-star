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

/*
void calc_density_for_all(Vec_1d<Particle>& particles, Vec_1d<Particle>* grid, double h, int N_grid_cells) {
	int N = particles.size();

	for (int i = 0; i < N; i++) {
		particles[i].change_density(particles[i].get_m() * particles[i].W(particles[i], h));
		int n_grid = N_grid_cells / 2 + particles[i].get_position()[0] / (2 * h);
		int m_grid =  N_grid_cells / 2 - particles[i].get_position()[1] / (2 * h);
		for (int k = m_grid - 1; k < m_grid + 2; k++) {
			if (k < 0 || k > N_grid_cells - 1) {
				continue;
			}
			for (int l = n_grid - 1; l < n_grid + 2; l++) {
				if (l < 0 || l > N_grid_cells - 1) {
					continue;
				}
				for (Particle particle : grid[N_grid_cells * k + l]) {
					double rho_ij = particles[i].get_m() * particles[i].W(particle, h);
					particles[i].change_density(particles[i].get_density() + rho_ij);
					particle.change_density(particle.get_density() + rho_ij);
				}
			}
		}
	}
}

void calc_pressure_for_all(Vec_1d<Particle>& particles, double k, double n) {
	int N = particles.size();
	for (int i = 0; i < N; i++) {
        particles[i].change_pressure(k * pow(particles[i].get_density(), 1 + 1 / n));
    }
}


Vec_2d<double> acceleration(Vec_1d<Particle> particles, Vec_1d<Particle>* grid, double h, double k, double n, double lmbda, double nu, int N_grid_cells) {
    int N = particles.size();
    Vec_2d<double> a(N, Vec_1d<double>(2));
    calc_density_for_all(particles, grid, h, N_grid_cells);
    calc_pressure_for_all(particles, k, n);
    Vec_1d<double> a_p(2);
    for (int i = 0; i < N; i++) {
        a[i][0] += - nu * particles[i].get_velocity()[0] - lmbda * particles[i].get_position()[0];
        a[i][1] += - nu * particles[i].get_velocity()[1] - lmbda * particles[i].get_position()[1];
        for (int j = i + 1; j < N; j++) {
            a_p[0] = - particles[i].get_m() * (particles[i].get_pressure() / (particles[i].get_density() * particles[i].get_density()) +
											 particles[j].get_pressure() / (particles[j].get_density() * particles[j].get_density())) * particles[i].gradW(particles[j], h)[0];
            a_p[1] = - particles[i].get_m() * (particles[i].get_pressure() / (particles[i].get_density() * particles[i].get_density()) +
						           particles[j].get_pressure() / (particles[j].get_density() * particles[j].get_density())) * particles[i].gradW(particles[j], h)[1];
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
    int Nt = tmax / dt;
    Vec_2d<Particle> result (Nt);

    const int N_grid_cells = d / (2 * h);
    Vec_1d<Particle> grid[N_grid_cells * N_grid_cells];

    for(int i = 0; i < N; i ++){
        result[0].push_back(particles[i]);
        int n_grid = N_grid_cells / 2 + particles[i].get_position()[0] / (2 * h);
        int m_grid =  N_grid_cells / 2 - particles[i].get_position()[1] / (2 * h);
        grid[N_grid_cells * m_grid + n_grid].push_back(particles[i]);
    }

    Vec_2d<double> acc = acceleration(particles, grid, h, k, n, lmbda, nu, N_grid_cells);
    Vec_1d<double> rho (N);
    for (int i = 0; i < Nt; i ++){
        for (int j = 0; j < N; j++){
            particles[j].change_velocity(acc[j], dt);
            particles[j].change_position(particles[j].get_velocity(), dt);
            acc = acceleration(particles, grid, h, k, n, lmbda, nu, N_grid_cells);
            particles[j].change_velocity(acc[j], dt);
        }

				for (int k = 0; k < N_grid_cells * N_grid_cells; k++) {
	          grid[k].clear();
	      }

        for (int j = 0; j < N; j ++) {
            result[i].push_back(particles[j]);
            int n_grid =  N_grid_cells / 2 + particles[j].get_position()[0] / (2 * h);
            int m_grid =  N_grid_cells / 2 - particles[j].get_position()[1] / (2 * h);
            grid[N_grid_cells * m_grid + n_grid].push_back(particles[j]);
        }
    }
    return result;
}
*/

int main() {
    const int N_part = 50;
    double h = 0.1;
    double d = 8.;
    Vec_2d<Particle> result;
    Vec_1d<Particle> particles;
    double k = 0.1, n = 1, nu = 1, lmbda = 4, dt = 0.04, maxt = 6;
    Vec_2d<double> a(N_part, Vec_1d<double>(2));

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> NormRand(0, 0.1);

    {
        Timer<std::chrono::seconds> t;

        for (int i = 0; i < N_part; i++) {
            particles.push_back(Particle(Vec_1d<double>{NormRand(gen), NormRand(gen), 0},
                                         Vec_1d<double>{NormRand(gen), NormRand(gen), 0},
                                         0.03)
            );
        }

        result = calc(particles, h, d, k, n, lmbda, nu, maxt, dt);
    }

    std::vector<sf::CircleShape> sprites(N_part);
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
