/**
 * @brief: generates a cuboid of particles specified by the user input
 */
#include "ParticleGenerator.h"

void ParticleGenerator::instantiateCuboid(ParticleContainer &container, std::array<double, 3> llfc,
                                          std::array<unsigned int, 3> particlePerDimension, std::array<double, 3> particleVelocity,
                                          double h, double mass, int generateNumber)
{
    double meanValueVelocity{0.1};
    std::array<double, 3> mbVelocity = MaxwellBoltzmannDistribution::maxwellBoltzmannDistributedVelocity(meanValueVelocity, 3); 
    auto amountOfParticles = particlePerDimension[0] * particlePerDimension[1] * particlePerDimension[2];
    std::vector<Particle> particles = container.getParticles();
    particles.reserve(particles.size() + amountOfParticles); // reserve space for amountOfParticles in the vector for performance
    for (unsigned int i = 0; i < particlePerDimension[0]; ++i)
    {
        for (unsigned int j = 0; j < particlePerDimension[1]; ++j)
        {
            for (unsigned int k = 0; k < particlePerDimension[2]; ++k)
            {
                std::array<double, 3> x_arg{i * h + llfc[0], j * h + llfc[1], k * h + llfc[2]};
                std::array<double, 3> v_arg{particleVelocity + mbVelocity}; // Calculate via the Brownian motion
                particles.emplace_back(Particle{x_arg, v_arg, mass, generateNumber});
            }
        }
    }
    container.setParticles(particles);
}

/* Use this constructor for the initialization
Particle::Particle(std::array<double, 3> x_arg, std::array<double, 3> v_arg, double m_arg, int type_arg):x(x_arg), v(v_arg), m(m_arg), type(type_arg) {
    f = {0., 0., 0.};
    old_f = {0., 0., 0.};
    std::cout << "Particle generated!" << std::endl;
}*/
