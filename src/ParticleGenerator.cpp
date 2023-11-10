/**
 * @brief: generates a cuboid of particles specified by the user input
 */
#include "ParticleGenerator.h"

void ParticleGenerator::instantiateCuboid(ParticleContainer &container, std::array<double, 3> llfc,
                                          std::array<unsigned int, 3> particlePerDimension, double h, double mass,
                                          std::array<double, 3> particleVelocity) {
    double meanValueVelocity{0.0}; // TODO: Is it the Erwartungswert of the normal distribution
    std::array<double, 3> mbVelocity = maxwellBoltzmannDistributedVelocity(meanValueVelocity, 3);
    std::vector<Particle> particles_temp; //TODO: reserve space for the vector for Performance
    for (int i = 0; i < particlePerDimension[0]; ++i) {
        for (int j = 0; j < particlePerDimension[1]; ++j) {
            for (int k = 0; k < particlePerDimension[2]; ++k) {
                std::array<double, 3>  x_arg{i*h, j*h, k*h};
                auto v_arg = particleVelocity + mbVelocity;
                //std::array<double, 3> v_arg{particleVelocity + mbVelocity}; // Calculate via the Brownian motion
                particles_temp.emplace_back(Particle{x_arg, v_arg, mass, 0});
            }
        }
    }
    container.setParticles(particles_temp);
}

/* Use this constructor for the initialization
Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg):x(x_arg), v(v_arg), m(m_arg), type(type_arg) {
    f = {0., 0., 0.};
    old_f = {0., 0., 0.};
    std::cout << "Particle generated!" << std::endl;
}*/
