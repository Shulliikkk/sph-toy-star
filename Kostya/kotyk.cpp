#include <SFML/Graphics.hpp>
#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <random>
#include <unistd.h>
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
	float rho;
	float p;
	Vec_1d<float> velocity =Vec_1d<float>(2);
	Vec_1d<float> position = Vec_1d<float>(2);

public:
  Particle(Vec_1d<float> position, Vec_1d<float> velocity, float m): m(m), velocity(velocity), position(position){
    }

  float distance(Particle& other){
        return std::sqrt(std::pow(this->position[0]-other.position[0], 2) +
                         std::pow(this->position[1]-other.position[1], 2));
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
                              scalar_deriv * (this->position[1] - other.position[1])};
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
        for (int k = 0; k < 2; k ++){
            velocity[k] += acc[k] * dt / 2;
        }
    }

    void change_position(Vec_1d<float> velocity, float dt){
        for (int k = 0; k < 2; k ++){
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

void calc_density_for_all(Vec_1d<Particle>& particles, Vec_1d<Particle>* grid, float h, int N_grid_cells) {
	int N = particles.size();

	for(int i = 0; i < N; i++) {
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
					float rho_ij = particles[i].get_m() * particles[i].W(particle, h);
					particles[i].change_density(particles[i].get_density() + rho_ij);
					particle.change_density(particle.get_density() + rho_ij);
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
    Vec_1d<float> rho(N);
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

        for (int j = 0; j < N; j++) {
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
		    const int N_part = 20;
		    float h = 0.1;
		    float d = 8;

		    std::vector <sf::CircleShape> sprites(N_part);
		    std::random_device rd;
		    std::mt19937 gen(rd());
		    std::normal_distribution<> NormRand(0, 0.1);

		    Vec_1d<Particle> particles;
		    for (int i = 0; i < N_part; i++) {
		        particles.push_back(Particle(Vec_1d<float>{NormRand(gen), NormRand(gen), 0}, Vec_1d<float>{NormRand(gen),  NormRand(gen), 0}, 0.03));
		    }

		    float k = 0.1, n = 1, nu = 1, lmbda = 4, dt = 0.04, maxt = 6;
		    Vec_2d<float> a(N_part, Vec_1d<float>(2));
		    Vec_2d<Particle> result = calc(particles, h, d, k, n, lmbda, nu, maxt, dt);
		    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
		    std::cout << "Diff(ms) = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << std::endl;

        std::size_t step = 0;
        sf::RenderWindow window (sf::VideoMode (400, 400), " SFML works!");
        float want_fps = 5;
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
            for(std::size_t j = 0; j < result[step].size(); j++){
              sprite.setPosition(result[step][j].get_position()[0] * 100 + 200, result[step][j].get_position()[1] * 60 + 200);
              window.draw(sprite);
            }
          }

            window.display();
            if (step + 1 < result.size())
                step += 1;
            else
                step = 0;
            sf::Int32 frame_duration = loop_timer.getElapsedTime().asMilliseconds();
        sf::Int32 time_to_sleep = int(1000.f/want_fps) - frame_duration;
        if (time_to_sleep > 0) {
            sf::sleep(sf::milliseconds(time_to_sleep));
        }
        loop_timer.restart();
        }
}
