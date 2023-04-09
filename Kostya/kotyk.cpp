#include <SFML/Graphics.hpp>
#include<vector>
#include<iostream>
#include<cmath>
#include<iomanip>
#include <random>
#include<unistd.h>
#include <chrono>

template<typename T>
using Vec_1d = std::vector<T>;
template<typename T>
using Vec_2d = std::vector<std::vector<T>>;
template<typename T>
using Vec_3d = std::vector<std::vector<std::vector<T>>>;

#define float double

const float PI = std::acos(-1.0);;

class Particle{
private:
	float m;
	Vec_1d<float> velocity =Vec_1d<float>(3);
	Vec_1d<float> position = Vec_1d<float>(3);
	float rho;
	float p;

public:
    Particle(Vec_1d<float> position, Vec_1d<float> velocity, float m): m(m), velocity(velocity), position(position){
    }

    float distance(Particle& other){
        return std::sqrt(std::pow(this->position[0]-other.position[0], 2) +
                         std::pow(this->position[1]-other.position[1], 2) +
                         std::pow(this->position[2]-other.position[2], 2));
    }
	float W(Particle& other, float h){
	    float abs_r = distance(other);
        return 1 / (h * h * h + std::pow(PI, 1.5))
                * std::exp(- abs_r * abs_r / (h * h));
	}

	Vec_1d<float> gradW(Particle& other, float h){
	    float abs_r = distance(other);
	    float scalar_deriv = -2 / (h * h * h * h * h * std::pow(PI, 1.5))
                                * std::exp(- abs_r * abs_r / (h * h));
	    Vec_1d<float> grad = {scalar_deriv * (this->position[0] - other.position[0]),
                              scalar_deriv * (this->position[1] - other.position[1]),
                              scalar_deriv * (this->position[2] - other.position[2])};
        return grad;
	}

    Vec_1d<float> get_position(){
        return position;
    }

    Vec_1d<float> get_velocity(){
        return velocity;
    }

    float get_m(){
        return m;
    }

		float get_density() {
			return rho;
		}

		float get_pressure() {
			return p;
		}

    void change_velocity(Vec_1d<float>& acc, float dt){
        for (int k = 0; k < 3; k ++){
            velocity[k] += acc[k] * dt / 2;
        }
    }

    void change_position(Vec_1d<float> velocity, float dt){
        for (int k = 0; k < 3; k ++){
            position[k] += velocity[k] * dt;
        }
    }

		void change_density(float density) {
			rho = density;
		}

		void change_pressure(float pressure) {
			p = pressure;
		}
};

/*
void calc_density_for_all(Vec_1d<Particle>& particles, Vec_1d<Particle>* grid, float h, int N_grid_cells) {
	int N = particles.size();

	for(int i = 0; i < N; i++) {
		particles[i].change_density(particles[i].get_m() * particles[i].W(particles[i], h));
		int n_grid = N_grid_cells / 2 + particles[i].get_position()[0] / (2 * h);
		int m_grid =  N_grid_cells / 2 - particles[i].get_position()[1] / (2 * h);
		//int n_grid = particles[i].get_position()[0] / (2 * h);
		//int m_grid = particles[i].get_position()[1] / (2 * h);
		for (int k = m_grid - 1; k < m_grid + 2; k++) {
			if (k < 0 || k > N_grid_cells - 1) {
				continue;
			}
			for (int l = n_grid - 1; l < n_grid + 2; l++) {
				if (l < 0 || l > N_grid_cells - 1) {
					continue;
				}
				for (Particle particle : grid[N_grid_cells * k + l]) {
					float rho_ij = particles[i].get_m() * particles[i].W(particle, h);
					particles[i].change_density(particles[i].get_density() + rho_ij);
					particle.change_density(particle.get_density() + rho_ij);
				}
			}
		}
	}
}
*/

/*
void calc_density_for_all(Vec_1d<Particle>& particles, Vec_1d<Particle>* grid, float h, int N_grid_cells) {
	int N = particles.size();
	for(int i = 0; i < N; i++) {
		particles[i].change_density(particles[i].get_m() * particles[i].W(particles[i], h));
		int n_grid = particles[i].get_position()[0] / (2 * h);
		int m_grid = particles[i].get_position()[1] / (2 * h);
		for (int k = m_grid - 1; k < m_grid + 2; k++) {
			if (k < 0 || k > N_grid_cells - 1) {
				continue;
			}
			for (int l = n_grid - 1; l < n_grid + 2; l++) {
				if (l < 0 || l > N_grid_cells - 1) {
					continue;
				}
				for (Particle particle : grid[N_grid_cells * k + l]) {
					float rho_ij = particles[i].get_m() * particles[i].W(particle, h);
					particles[i].change_density(particles[i].get_density() + rho_ij);
					particle.change_density(particle.get_density() + rho_ij);
				}
			}
		}
	}
}
*/


void calc_density_for_all(Vec_1d<Particle>& particles, Vec_1d<Particle>* grid, float h, int N_grid_cells) {
	//std::cout << N_grid_cells << '\n';
	for (int i = 0; i < N_grid_cells; i++) {
		for (int j = 0; j < N_grid_cells; j++) {
			//std::cout << "kek " << grid[N_grid_cells * i + j].size() << '\n';
			for (Particle particle_1 : grid[N_grid_cells * i + j]) {
				for (Particle particle_2 : grid[N_grid_cells * i + j]) {
					float rho_ij = particle_1.get_m() * particle_1.W(particle_2, h);
					particle_1.change_density(particle_1.get_density() + rho_ij);
				}
				//std::cout << "lol " << i << ' ' << j << '\n';
				if (j + 1 < N_grid_cells) {
					for (Particle particle_2 : grid[N_grid_cells * i + (j + 1)]) {
						float rho_ij = particle_1.get_m() * particle_1.W(particle_2, h);
						particle_1.change_density(particle_1.get_density() + rho_ij);
						particle_2.change_density(particle_2.get_density() + rho_ij);
					}
					//std::cout << "j + 1 " << i << ' ' << j << '\n';
				}

				if (j + 1 < N_grid_cells && 0 < i - 1) {
					for (Particle particle_2 : grid[N_grid_cells * (i - 1) + (j + 1)]) {
						float rho_ij = particle_1.get_m() * particle_1.W(particle_2, h);
						particle_1.change_density(particle_1.get_density() + rho_ij);
						particle_2.change_density(particle_2.get_density() + rho_ij);
					}
					//std::cout << "i + 1 " << i << ' ' << j << '\n';
				}

				if (i + 1 < N_grid_cells) {
					for (Particle particle_2 : grid[N_grid_cells * (i + 1) + j]) {
						float rho_ij = particle_1.get_m() * particle_1.W(particle_2, h);
						particle_1.change_density(particle_1.get_density() + rho_ij);
						particle_2.change_density(particle_2.get_density() + rho_ij);
					}
					//std::cout << "i + 1 " << i << ' ' << j << '\n';
				}

				if ((j + 1 < N_grid_cells) && (i + 1 < N_grid_cells)) {
					for (Particle particle_2 : grid[N_grid_cells * (i + 1) + (j + 1)]) {
						float rho_ij = particle_1.get_m() * particle_1.W(particle_2, h);
						particle_1.change_density(particle_1.get_density() + rho_ij);
						particle_2.change_density(particle_2.get_density() + rho_ij);
					}
					//std::cout << "j + 1; i + 1 " << i << ' ' << j << '\n';
				}
			}
		}
	}
}


void calc_pressure_for_all(Vec_1d<Particle>& particles, float k, float n) {
	int N = particles.size();
	for (int i = 0; i < N; i++) {
    particles[i].change_pressure(k * pow(particles[i].get_density(), 1 + 1 / n));
  }
}


Vec_2d<float> acceleration(Vec_1d<Particle> particles, Vec_1d<Particle>* grid, float h, float k, float n, float lmbda, float nu, int N_grid_cells) {
    int N = particles.size();
    Vec_2d<float> a(N, Vec_1d<float>(2));
    calc_density_for_all(particles, grid, h, N_grid_cells);
    calc_pressure_for_all(particles, k, n);
    Vec_1d<float> a_p(2);
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

Vec_2d<Particle> calc(Vec_1d<Particle>& particles, float h, float d, float k, float n, float lmbda, float nu, float tmax, float dt){
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
    Vec_2d<float> acc = acceleration(particles, grid, h, k, n, lmbda, nu, N_grid_cells);
    Vec_1d<float> rho (N);
    for (int i = 0; i < Nt; i ++){
        for (int j = 0; j < N; j++){
            particles[j].change_velocity(acc[j], dt);
            particles[j].change_position(particles[j].get_velocity(), dt);
            acc = acceleration(particles, grid, h, k, n, lmbda, nu, N_grid_cells);
            particles[j].change_velocity(acc[j], dt);
        }

				/*
				for (int m = 0; m < N_grid_cells; m++) {
					for (int n = 0; n < N_grid_cells; n++) {
						Vec_1d<Particle> cell = grid[N_grid_cells * m + n];
						for (auto it_particle = cell.begin(); it_particle != cell.end(); it_particle++) {
							result[i].push_back(*it_particle);
							int n_grid = (*it_particle).get_position()[0] / (2 * h);
							int m_grid = (*it_particle).get_position()[1] / (2 * h);
							if (n_grid != n || m_grid != m) {
								cell.erase(it_particle);
								grid[N_grid_cells * m_grid + n_grid].push_back(*it_particle);
							}
						}
					}
				}
				*/


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

int main(){
		std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    float h = 0.1;
		float d = 6;
		Vec_1d<Particle> particles = {Particle(Vec_1d<float>{-1.54909388, 1.40804976, 0}, Vec_1d<float>{0.025548152,  0.082706075, 0}, 0.2),
																		Particle(Vec_1d<float>{1.43990639,  -1.59498159, 0}, Vec_1d<float>{0.04342804,   0.048848213,0}, 0.2),
																		Particle(Vec_1d<float>{-1.27890266, 1.26921543, 0}, Vec_1d<float>{0.081535258,  0.019112969, 0}, 0.2),
																		Particle(Vec_1d<float>{1.60406402,  1.32339721, 0}, Vec_1d<float>{0.094470003, 0.088040316, 0}, 0.2),
																		Particle(Vec_1d<float>{1.07615878, -1.62474234, 0}, Vec_1d<float>{0.049821641, 0.07471592, 0}, 0.2),
																		Particle(Vec_1d<float>{1.32821573,  1.00483002, 0}, Vec_1d<float>{0.056551748,  0.054349262, 0}, 0.2),
																		Particle(Vec_1d<float>{-1.34537595,  1.7395526 , 0}, Vec_1d<float>{0.006876223, 0.062816352, 0}, 0.2),
																		Particle(Vec_1d<float>{1.02683444,  1.73286105, 0}, Vec_1d<float>{0.027522109, 0.035860894, 0}, 0.2),
																		Particle(Vec_1d<float>{-1.07508888,  1.96809424, 0}, Vec_1d<float>{0.018117553,  0.016827289, 0}, 0.2),
																		Particle(Vec_1d<float>{1.42200362, -1.61358428, 0}, Vec_1d<float>{0.075123249,  0.06037533, 0}, 0.2),
																		Particle(Vec_1d<float>{1.94909388, 1.20804976, 0}, Vec_1d<float>{0.025548152,  0.082706075, 0}, 0.2)};

    float k = 0.1, n = 1, nu = 1, lmbda = 4, dt = 0.0005, maxt = 0.4;
    Vec_2d<float> a(10, Vec_1d<float>(2));
    Vec_2d<Particle> result = calc(particles, h, d, k, n, lmbda, nu, maxt, dt);
    std::cout << std::setprecision(7);

    for (std::size_t i = 0; i < result.size(); i++){
        for(std::size_t j = 0; j < result[i].size(); j++){
            std::cout << result[i][j].get_position()[0] << ' ' << result[i][j].get_position()[1] << '\n';
        }
    }
		std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		std::cout << "Diff(ms) = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << std::endl;
}
