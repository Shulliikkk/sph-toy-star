#include <SFML/Graphics.hpp>
#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <random>
#include <unistd.h>
#include <algorithm>

template<typename T>
using Vec_1d = std::vector<T>;

template<typename T>
using Vec_2d= std::vector<std::vector<T>>;

template<typename T>
using Vec_3d = std::vector<std::vector<std::vector<T>>>;

#define float double

const float PI = std::acos(-1.0);

class Particle {
private:
	float m;
	float rho;
	float p;
	Vec_1d<float> velocity =Vec_1d<float>(2);
	Vec_1d<float> position = Vec_1d<float>(2);

public:

  Particle(Vec_1d<float> position, Vec_1d<float> velocity, float m): m(m), velocity(velocity), position(position) {}

  float distance(Particle& other) {
        return std::sqrt(std::pow(this->position[0]-other.position[0], 2) +
                         std::pow(this->position[1]-other.position[1], 2));
  }

	float W(Particle& other, float h) {
	    float abs_r = distance(other);
        return 1 / (h * h * h + std::pow(PI, 1.5))
                * std::exp(- abs_r * abs_r / (h * h));
	}

	Vec_1d<float> gradW(Particle& other, float h) {
	    float abs_r = distance(other);
	    float scalar_deriv = -2 / (h * h * h * h * h * std::pow(PI, 1.5))
                                * std::exp(- abs_r * abs_r / (h * h));
	    Vec_1d<float> grad = {scalar_deriv * (this -> position[0] - other.position[0]),
                              scalar_deriv * (this -> position[1] - other.position[1])};
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

class Processing {
private:
  Vec_2d<float> density;

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
      Vec_1d<float> curr_density(N);
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
          curr_density[i] = particles[i].get_density();
      }
      density.push_back(curr_density);
      return a;
  }

public:
  Processing() {};

  Vec_2d<Particle> calc(int N, float h, float d, float k, float n, float lmbda, float nu, float tmax, float dt){
      std::random_device rd;
      std::mt19937 gen(rd());
      std::normal_distribution<> NormRand(0, 0.2);

      Vec_1d<Particle> particles;

      for (int i = 0; i < N; i++) {
          particles.push_back(Particle(Vec_1d<float>{NormRand(gen), NormRand(gen), 0}, Vec_1d<float>{NormRand(gen),  NormRand(gen), 0}, 0.03));
      }

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
      for (int i = 1; i < Nt; i ++){
          for (int j = 0; j < N; j++){
              particles[j].change_velocity(acc[j], dt);
              particles[j].change_position(particles[j].get_velocity(), dt);
          }
          acc = acceleration(particles, grid, h, k, n, lmbda, nu, N_grid_cells);

  				for (int k = 0; k < N_grid_cells * N_grid_cells; k++) {
  	          grid[k].clear();
  	      }

          for (int j = 0; j < N; j++) {
              result[i].push_back(particles[j]);
              int n_grid =  N_grid_cells / 2 + particles[j].get_position()[0] / (2 * h);
              int m_grid =  N_grid_cells / 2 - particles[j].get_position()[1] / (2 * h);
              grid[N_grid_cells * m_grid + n_grid].push_back(particles[j]);}
      }
      return result;
  	}

    Vec_2d<float> get_density_for_all() {
      return density;
    }
};

class Visualisation {
private:
  sf::Color hsv(int hue, float sat, float val) {
      hue %= 360;
      while(hue < 0) hue += 360;

      if(sat < 0.f) sat = 0.f;
      if(sat > 1.f) sat = 1.f;

      if(val < 0.f) val = 0.f;
      if(val > 1.f) val = 1.f;

      int h = hue / 60;
      float f = float(hue) / 60 - h;
      float p = val * (1.f - sat);
      float q = val * (1.f - sat * f);
      float t = val * (1.f - sat * (1 - f));

      switch(h) {
          default:
          case 0:
          case 6: return sf::Color(val * 255, t * 255, p * 255);
          case 1: return sf::Color(q * 255, val * 255, p * 255);
          case 2: return sf::Color(p * 255, val * 255, t * 255);
          case 3: return sf::Color(p * 255, q * 255, val * 255);
          case 4: return sf::Color(t * 255, p * 255, val * 255);
          case 5: return sf::Color(val * 255, p * 255, q * 255);
      }
  }

public:
  Visualisation() {};
  void loop(Vec_2d<Particle> result, Vec_2d<float> density, float h) {
    int N_part = result[0].size();
    std::vector <sf::CircleShape> sprites(N_part);
    std::vector <float> rhos;
    rhos.clear();

    for (std::size_t i = 0; i < result.size(); i++) {
        for(std::size_t j = 0; j < result[i].size(); j++) {
            float rho = density[i][j];
            rhos.push_back(rho);
        }
    }

    sf::RenderWindow window (sf::VideoMode (600, 600), "Toy Star Simulation");
    float want_fps = 60;
    sf::Clock loop_timer;

    auto it = std::minmax_element(rhos.begin(), rhos.end(), std::greater<float>());
    auto max_rho = *it.first;
    auto min_rho = *it.second;

    float rhos_range = max_rho - min_rho;
    float rhos_step = rhos_range / 360;

    float rhos_quarter = min_rho + 360 * rhos_step / 4;
    float rhos_half = rhos_quarter + 180 * rhos_step / 2;

    std::sort(rhos.begin(), rhos.end() - 1);

    std::size_t step = 0;
    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear();

        for (int i = 0; i < sprites.size(); i++) {
            sprites[i].setRadius(2);
            sprites[i].setOutlineThickness(1.5);

            for(std::size_t j = 0; j < result[step].size(); j++) {
                sprites[i].setPosition(result[step][j].get_position()[0] * 100 + 300, result[step][j].get_position()[1] * 100 + 300);
                float rho = density[step][j];
                if (rho < rhos_quarter)
                {
                    for (int k = 1; k < 360; k++)
                    {
                        if (min_rho + (k - 1) * rhos_step / 4 <= rho && rho < min_rho + k * rhos_step / 4)
                        {
                            sprites[i].setFillColor(hsv(k, 1.f, 0.5));
                            sprites[i].setOutlineColor(hsv(k, 1.f, 1.f));
                        }
                    }
                }

                else if (rho < rhos_half)
                {
                    for (int k = 1; k < 180; k++)
                    {
                        if (rhos_quarter + (k - 1) * rhos_step / 2 <= rho && rho < rhos_quarter + k * rhos_step / 2)
                        {
                            sprites[i].setFillColor(hsv(k, 1.f, 0.5));
                            sprites[i].setOutlineColor(hsv(k, 1.f, 1.f));
                        }
                    }
                }

                else
                {
                    for (int k = 1; k < 180; k++)
                    {
                        if (rhos_half + (k - 1) * rhos_step <= rho && rho < rhos_half + k * rhos_step)
                        {
                            sprites[i].setFillColor(hsv(k, 1.f, 0.5));
                            sprites[i].setOutlineColor(hsv(k, 1.f, 1.f));
                        }
                    }
                }

                window.draw(sprites[i]);
            }
        }


        //шкала
        sf::VertexArray quad1(sf::Quads, 4);
        sf::VertexArray quad2(sf::Quads, 4);
        sf::VertexArray quad3(sf::Quads, 4);

        quad1[0].position = sf::Vector2f(560.f, 155.f);
        quad1[1].position = sf::Vector2f(590.f, 155.f);
        quad1[2].position = sf::Vector2f(590.f, 255.f);
        quad1[3].position = sf::Vector2f(560.f, 255.f);

        quad2[0].position = sf::Vector2f(560.f, 255.f);
        quad2[1].position = sf::Vector2f(590.f, 255.f);
        quad2[2].position = sf::Vector2f(590.f, 350.f);
        quad2[3].position = sf::Vector2f(560.f, 350.f);

        quad3[0].position = sf::Vector2f(560.f, 350.f);
        quad3[1].position = sf::Vector2f(590.f, 350.f);
        quad3[2].position = sf::Vector2f(590.f, 445.f);
        quad3[3].position = sf::Vector2f(560.f, 445.f);


        quad1[0].color = sf::Color{255, 0, 255, 255};
        quad1[1].color = sf::Color{255, 0, 255, 255};
        quad1[2].color = sf::Color{0, 255, 255, 255};
        quad1[3].color = sf::Color{0, 255, 255, 255};

        quad2[0].color = sf::Color{0, 255, 255, 255};
        quad2[1].color = sf::Color{0, 255, 255, 255};
        quad2[2].color = sf::Color{255, 255, 0, 255};
        quad2[3].color = sf::Color{255, 255, 0, 255};

        quad3[0].color = sf::Color{255, 255, 0, 255};
        quad3[1].color = sf::Color{255, 255, 0, 255};
        quad3[2].color = sf::Color{255, 0, 0, 255};
        quad3[3].color = sf::Color{255, 0, 0, 255};

       sf::Text text_min;
       sf::Text text_max;
       sf::Text text_half;
       sf::Text text_quarter;

       sf::Font font;
       font.loadFromFile("arial.ttf");

       text_min.setFont(font);
       text_max.setFont(font);
       text_half.setFont(font);
       text_quarter.setFont(font);

       std::string min_str = std::to_string(min_rho);
       std::string max_str = std::to_string(max_rho);
       std::string half_str = std::to_string(rhos_half);
       std::string quarter_str = std::to_string(rhos_quarter);

       text_min.setString(min_str);
       text_max.setString(max_str);
       text_half.setString(half_str);
       text_quarter.setString(quarter_str);

       int size = 12;

       text_min.setCharacterSize(size);
       text_max.setCharacterSize(size);
       text_half.setCharacterSize(size);
       text_quarter.setCharacterSize(size);

       text_min.setFillColor(sf::Color::White);
       text_max.setFillColor(sf::Color::White);
       text_half.setFillColor(sf::Color::White);
       text_quarter.setFillColor(sf::Color::White);

       text_min.setPosition(560.f - 8 * min_str.length(), 442.f - size / 2);
       text_max.setPosition(560.f - 8 * max_str.length(), 152.f - size / 2);
       text_half.setPosition(560.f - 8 * half_str.length(), 297.f - size / 2);
       text_quarter.setPosition(560.f - 8 * quarter_str.length(), 369.5 - size / 2);

       sf::Vertex line1[] =
       {
           sf::Vertex(sf::Vector2f(560.f, 155.f)),
           sf::Vertex(sf::Vector2f(590.f, 155.f))
       };

       sf::Vertex line2[] =
       {
           sf::Vertex(sf::Vector2f(590.f, 155.f)),
           sf::Vertex(sf::Vector2f(590.f, 445.f))
       };

       sf::Vertex line3[] =
       {
           sf::Vertex(sf::Vector2f(590.f, 445.f)),
           sf::Vertex(sf::Vector2f(560.f, 445.f))
       };

       sf::Vertex line4[] =
       {
           sf::Vertex(sf::Vector2f(560.f, 445.f)),
           sf::Vertex(sf::Vector2f(560.f, 155.f))
       };

       sf::Vertex line5[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 300.f)),
           sf::Vertex(sf::Vector2f(590.f, 300.f))
       };

       sf::Vertex line6[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 372.5)),
           sf::Vertex(sf::Vector2f(590.f, 372.5))
       };

       sf::Vertex line7[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 155.f)),
           sf::Vertex(sf::Vector2f(590.f, 155.f))
       };

       sf::Vertex line8[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 445.f)),
           sf::Vertex(sf::Vector2f(590.f, 445.f))
       };

       sf::Vertex line9[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 299.f)),
           sf::Vertex(sf::Vector2f(590.f, 299.f))
       };

       sf::Vertex line10[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 371.5)),
           sf::Vertex(sf::Vector2f(590.f, 371.5))
       };

       sf::Vertex line11[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 154.f)),
           sf::Vertex(sf::Vector2f(590.f, 154.f))
       };

       sf::Vertex line12[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 444.f)),
           sf::Vertex(sf::Vector2f(590.f, 444.f))
       };
       sf::Vertex line13[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 301.f)),
           sf::Vertex(sf::Vector2f(590.f, 301.f))
       };

       sf::Vertex line14[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 373.5)),
           sf::Vertex(sf::Vector2f(590.f, 373.5))
       };

       sf::Vertex line15[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 156.f)),
           sf::Vertex(sf::Vector2f(590.f, 156.f))
       };

       sf::Vertex line16[] =
       {
           sf::Vertex(sf::Vector2f(558.f, 446.f)),
           sf::Vertex(sf::Vector2f(590.f, 446.f))
       };

       float graphic_step = 290 / 360;

       sf::Vertex line17[] =
       {
           sf::Vertex(sf::Vector2f(560.f, 445.f - graphic_step / 4)),
           sf::Vertex(sf::Vector2f(590.f, 445.f - graphic_step / 4))
       };

       sf::Vertex line18[] =
       {
           sf::Vertex(sf::Vector2f(560.f, 445.f - graphic_step / 2)),
           sf::Vertex(sf::Vector2f(590.f, 445.f - graphic_step / 2))
       };
       sf::Vertex line19[] =
       {
           sf::Vertex(sf::Vector2f(560.f, 372.5 - graphic_step / 2)),
           sf::Vertex(sf::Vector2f(590.f, 372.5 - graphic_step / 2))
       };

       sf::Vertex line20[] =
       {
           sf::Vertex(sf::Vector2f(560.f, 372.5 - graphic_step)),
           sf::Vertex(sf::Vector2f(590.f, 372.5 - graphic_step))
       };

       sf::Vertex line21[] =
       {
           sf::Vertex(sf::Vector2f(560.f, 300.f - graphic_step)),
           sf::Vertex(sf::Vector2f(590.f, 300.f - graphic_step))
       };

       sf::Vertex line22[] =
       {
           sf::Vertex(sf::Vector2f(560.f, 300.f - 2 * graphic_step)),
           sf::Vertex(sf::Vector2f(590.f, 300.f - 2 * graphic_step))
       };

       window.draw(text_min);
       window.draw(text_max);
       window.draw(text_half);
       window.draw(text_quarter);

       window.draw(quad1);
       window.draw(quad2);
       window.draw(quad3);
       window.draw(line1, 2, sf::Lines);
       window.draw(line2, 2, sf::Lines);
       window.draw(line3, 2, sf::Lines);
       window.draw(line4, 2, sf::Lines);
       window.draw(line5, 2, sf::Lines);
       window.draw(line6, 2, sf::Lines);
       window.draw(line7, 2, sf::Lines);
       window.draw(line8, 2, sf::Lines);
       window.draw(line9, 2, sf::Lines);
       window.draw(line10, 2, sf::Lines);
       window.draw(line11, 2, sf::Lines);
       window.draw(line12, 2, sf::Lines);
       window.draw(line13, 2, sf::Lines);
       window.draw(line14, 2, sf::Lines);
       window.draw(line15, 2, sf::Lines);
       window.draw(line16, 2, sf::Lines);
       window.draw(line17, 2, sf::Lines);
       window.draw(line18, 2, sf::Lines);
       window.draw(line19, 2, sf::Lines);
       window.draw(line20, 2, sf::Lines);
       window.draw(line21, 2, sf::Lines);
       window.draw(line22, 2, sf::Lines);

       window.display();

        if (step + 1  < result.size())
            step += 1;
        else
            step = 0;

        sf::Int32 frame_duration = loop_timer.getElapsedTime().asMilliseconds();
        sf::Int32 time_to_sleep = int(1000.f/want_fps) - frame_duration;

        if (time_to_sleep > 0)
        {
            sf::sleep(sf::milliseconds(time_to_sleep));
        }

        loop_timer.restart();
    }
  }
};

int main()
{
    const int N_part = 200;
    float h = 0.1, d = 20, k = 0.1, n = 1, nu = 1, lmbda = 4, dt = 0.04, maxt = 12;
    Processing processing;
    Visualisation visualisation;
    Vec_2d<Particle> result =  processing.calc(N_part, h, d, k, n, lmbda, nu, maxt, dt);
    Vec_2d<float> density = processing.get_density_for_all();
    visualisation.loop(result, density, h);
}
