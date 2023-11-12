/**
 * @brief: generates a cuboid of particles specified by the user input
 */
#include "ParticleGenerator.h"

void ParticleGenerator::instantiateCuboid(ParticleContainer &container, std::array<double, 3> llfc,
                                          std::array<unsigned int, 3> particlePerDimension, double h, double mass,
                                          std::array<double, 3> particleVelocity, int generateNumber) {
    double meanValueVelocity{0.0}; // TODO: Is it the Erwartungswert of the normal distribution
    std::array<double, 3> mbVelocity = MaxwellBoltzmannDistribution::maxwellBoltzmannDistributedVelocity(meanValueVelocity, 3); //TODO Fehlermeldung multiple def of MB
    auto amountOfParticles = particlePerDimension[0] * particlePerDimension[1] + particlePerDimension[2];
    std::vector<Particle> particles_temp;
    particles_temp.reserve(amountOfParticles); //reserve space for amountOfParticles for the vector for performance
    for (unsigned int i = 0; i < particlePerDimension[0]; ++i) {
        for (unsigned int j = 0; j < particlePerDimension[1]; ++j) {
            for (unsigned int k = 0; k < particlePerDimension[2]; ++k) {
                std::array<double, 3>  x_arg{i*h + llfc[0], j*h+llfc[1], k*h+llfc[3]};
                std::array<double, 3> v_arg{particleVelocity + mbVelocity}; // Calculate via the Brownian motion
                particles_temp.emplace_back(Particle{x_arg, v_arg, mass, generateNumber});
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
