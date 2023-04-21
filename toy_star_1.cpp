#include<vector>
#include<iostream>
#include<cmath>
#include<iomanip>
template<typename T>
using Vec_1d = std::vector<T>;
template<typename T>
using Vec_2d = std::vector<std::vector<T>>;
template<typename T>
using Vec_3d = std::vector<std::vector<std::vector<T>>>;

#define float double

const float PI = std::acos(-1.0);;

class Particle{
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
        return 1 / (std::pow(h, 3) + std::pow(PI,1.5))
                * std::exp(-std::pow(abs_r,2)/std::pow(h,2));
	}

	Vec_1d<float> gradW(Particle& other, float h){
	    float abs_r = distance(other);
	    float scalar_deriv = -2 / (std::pow(h, 5) * std::pow(PI, 1.5))
                                * std::exp(-std::pow(abs_r, 2) / std::pow(h, 2));
	    Vec_1d<float> grad = {scalar_deriv * (this->position[0] - other.position[0]),
                              scalar_deriv * (this->position[1] - other.position[1]),
                              scalar_deriv * (this->position[2] - other.position[2])};
        return grad;
	}
    float density(Vec_1d<Particle>& particles, float h) {
        int N = particles.size();
        Vec_1d<float> R_ij(2);
        float rho = 0;
        for (int j = 0; j < N; j++) {
            float rho_ij = particles[j].m * this->W(particles[j], h);
            rho += rho_ij;
        }
        return rho;
    }

    float pressure(float density, float k, float n) {
        float p;
        p = k * std::pow(density, 1 + 1 / n);
        return p;
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

private:
	float m;
	Vec_1d<float> velocity =Vec_1d<float>(3);
	Vec_1d<float> position = Vec_1d<float>(3);
};

Vec_2d<float> acceleration(Vec_1d<Particle> particles, float h, float k, float n, float lmbda, float nu) {
    int N = particles.size();
    Vec_2d<float> a(N, Vec_1d<float>(2));
    Vec_1d<float> rho (N);
    Vec_1d<float> p (N);
    for (int i = 0; i < N; i++){
        rho[i] = particles[i].density(particles, h);
        p[i] = particles[i].pressure(rho[i], k, n);
    }
    Vec_1d<float> R_ij(2);
    Vec_1d<float> a_p(2);
    for (int i = 0; i < N; i++) {
        a[i][0] += - nu * particles[i].get_velocity()[0] - lmbda * particles[i].get_position()[0];
        a[i][1] += - nu * particles[i].get_velocity()[1] - lmbda * particles[i].get_position()[1];
        for (int j = i + 1; j < N; j++) {
            a_p[0] = - particles[i].get_m() * (p[i] / pow(rho[i], 2) + p[j] / pow(rho[j], 2)) * particles[i].gradW(particles[j], h)[0];
            a_p[1] = - particles[i].get_m() * (p[i] / pow(rho[i], 2) + p[j] / pow(rho[j], 2)) * particles[i].gradW(particles[j], h)[1];
            a[i][0] += a_p[0];
            a[i][1] += a_p[1];
            a[j][0] += -a_p[0];
            a[j][1] += -a_p[1];
            }
    }
    return a;
}

Vec_2d<Particle> calc(Vec_1d<Particle>& particles, float h, float k, float n, float lmbda, float nu, float tmax, float dt){
    int N = particles.size();
    int Nt = tmax / dt;
    Vec_2d<Particle> result (Nt);
    for(int i = 0; i < N; i ++){
        result[0].push_back(particles[i]);
    }
    Vec_2d<float> acc = acceleration(particles, h, k, n, lmbda, nu);
    Vec_1d<float> rho (N);
    for (int i = 0; i < Nt; i ++){
        for (int j = 0; j < N; j++){
            particles[j].change_velocity(acc[j], dt);
            particles[j].change_position(particles[j].get_velocity(), dt);
            acc = acceleration(particles, h, k, n, lmbda, nu);
            particles[j].change_velocity(acc[j], dt);
        }
        for(int j = 0; j < N; j ++){
            result[i].push_back(particles[j]);
        }
    }
    return result;
}


/*
Vec2d_f acceleration(Vec_1d<Particle>& particles, float h, float k, float n, float lmbda, float nu) {
  float N = pos.size();
  Vec2d_f a(N, Vec1d_f(2));
  Vec1d_f rho = density(pos, m, h);
  Vec1d_f p = pressure(rho, k, n);
  Vec1d_f R_ij(2);
  Vec1d_f a_p(2);
  for (int i = 0; i < N; i++) {
    a[i][0] += - nu * vel[i][0] - lmbda * pos[i][0];
    a[i][1] += - nu * vel[i][1] - lmbda * pos[i][1];
    for (int j = i + 1; j < N; j++) {
      R_ij[0] = pos[i][0] - pos[j][0];
      R_ij[1] = pos[i][1] - pos[j][1];
      a_p[0] = - m * (p[i] / pow(rho[i], 2) + p[j] / pow(rho[j], 2)) * gradW(R_ij, h)[0];
      a_p[1] = - m * (p[i] / pow(rho[i], 2) + p[j] / pow(rho[j], 2)) * gradW(R_ij, h)[1];
      a[i][0] += a_p[0];
      a[i][1] += a_p[1];
      a[j][0] += -a_p[0];
      a[j][1] += -a_p[1];
    }
  }
  return a;
}
*/
int main(){
    float h = 0.1;
    Vec_1d<float> p1 = {1.3, 1.3, 0};
    Vec_1d<float> p2 = {1.2, 1.2, 0};
    //Particle pp1 {1, p1, p1};
    //Particle pp2 {1, p2, p2};
    Vec_1d<Particle> particles = {Particle(Vec_1d<float>{-1.54909388, -1.40804976, 0}, Vec_1d<float>{-1.25548152,  1.82706075, 0}, 0.2),
                                    Particle(Vec_1d<float>{-0.43990639,  0.59498159, 0}, Vec_1d<float>{0.4342804,   0.48848213,0}, 0.2),
                                    Particle(Vec_1d<float>{-0.27890266, -1.26921543, 0}, Vec_1d<float>{0.81535258,  0.19112969, 0}, 0.2),
                                    Particle(Vec_1d<float>{-0.60406402,  1.32339721, 0}, Vec_1d<float>{0.94470003, -0.88040316, 0}, 0.2),
                                    Particle(Vec_1d<float>{ 2.07615878, -0.62474234, 0}, Vec_1d<float>{-2.49821641, -1.7471592, 0}, 0.2),
                                    Particle(Vec_1d<float>{ 0.32821573,  0.00483002, 0}, Vec_1d<float>{0.56551748,  1.54349262, 0}, 0.2),
                                    Particle(Vec_1d<float>{-0.34537595,  0.7395526 , 0}, Vec_1d<float>{-1.06876223, -0.62816352, 0}, 0.2),
                                    Particle(Vec_1d<float>{-1.02683444,  0.73286105, 0}, Vec_1d<float>{1.27522109, -1.35860894, 0}, 0.2),
                                    Particle(Vec_1d<float>{-0.07508888,  0.96809424, 0}, Vec_1d<float>{-0.18117553,  0.16827289, 0}, 0.2),
                                    Particle(Vec_1d<float>{-0.42200362, -1.61358428, 0}, Vec_1d<float>{0.75123249,  1.6037533, 0}, 0.2)};

    float k = 0.1, n = 1, nu = 1, lmbda = 4, dt = 0.04, maxt = 1;
    Vec_2d<float> a(10, Vec_1d<float>(2));
//    rho = density(particles, h);
    Vec_2d<Particle> result = calc(particles, h, k,n, lmbda, nu, maxt, dt);
    std::cout << std::setprecision(7);

    for (std::size_t i = 0; i < result.size(); i++){
        for(std::size_t j = 0; j < result[i].size(); j++){
            std::cout << result[i][j].get_position()[0] << ' ' << result[i][j].get_position()[1] << '\n';
        }
    }


   // Vec_1d<float> grad = pp1.gradW(pp2, h);
    //std::cout << grad[0] << ' ' << grad[1] << ' ' << grad[2] << '\n';

}
