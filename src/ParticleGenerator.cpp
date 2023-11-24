/**
 * @brief: generates a cuboid of particles specified by the user input
 */
#include "ParticleGenerator.h"
/*
void ParticleGenerator::instantiateCuboid(ParticleContainer &container, std::array<double, 3> llfc,
                                          std::array<unsigned int, 3> particlePerDimension, std::array<double, 3> particleVelocity,
                                          double h, double mass, int generateNumber)
{
    double meanValueVelocity{0.1}; // TODO: This can als be passed as a parameter
    auto amountOfParticles = particlePerDimension[0] * particlePerDimension[1] * particlePerDimension[2];
    std::vector<Particle> particles = container.getParticles();
    particles.reserve(particles.size() + amountOfParticles); // reserve space for amountOfParticles in the vector for performance
    for (unsigned int i = 0; i < particlePerDimension[0]; ++i)
    {
        for (unsigned int j = 0; j < particlePerDimension[1]; ++j)
        {
            for (unsigned int k = 0; k < particlePerDimension[2]; ++k)
            {
                // Change: I think the mbVelocity has to be calculated for every particle
                std::array<double, 3> mbVelocity = MaxwellBoltzmannDistribution::maxwellBoltzmannDistributedVelocity(meanValueVelocity, 2); // TODO: Does here the mean of the Brownian Motion mean the same as the average velocity
                std::array<double, 3> x_arg{i * h + llfc[0], j * h + llfc[1], k * h + llfc[2]};
                std::array<double, 3> v_arg{particleVelocity + mbVelocity}; // Calculate via the Brownian motion
                particles.emplace_back(x_arg, v_arg, mass, generateNumber);
            }
        }
    }
    container.setParticles(particles);
}*/

void ParticleGenerator::instantiateCuboid(ParticleContainer &container, std::array<double, 3> llfc,
                                          std::array<unsigned int, 3> particlePerDimension, std::array<double, 3> particleVelocity,
                                          double h, double mass, int generateNumber)
{
    double meanValueVelocity{0.1}; // TODO: This can als be passed as a parameter
    auto amountOfParticles = particlePerDimension[0] * particlePerDimension[1] * particlePerDimension[2];
    //std::vector<Particle> particles = container.getParticles();
    //particles.reserve(particles.size() + amountOfParticles); // reserve space for amountOfParticles in the vector for performance
    //container.
    for (unsigned int i = 0; i < particlePerDimension[0]; ++i)
    {
        for (unsigned int j = 0; j < particlePerDimension[1]; ++j)
        {
            for (unsigned int k = 0; k < particlePerDimension[2]; ++k)
            {
                // Change: I think the mbVelocity has to be calculated for every particle
                std::array<double, 3> mbVelocity = MaxwellBoltzmannDistribution::maxwellBoltzmannDistributedVelocity(meanValueVelocity, 2); // TODO: Does here the mean of the Brownian Motion mean the same as the average velocity
                std::array<double, 3> x_arg{i * h + llfc[0], j * h + llfc[1], k * h + llfc[2]};
                std::array<double, 3> v_arg{particleVelocity + mbVelocity}; // Calculate via the Brownian motion
                Particle temp = Particle{x_arg, v_arg, mass, generateNumber};
                container.add(temp); //using the add function in order to be able to add elements to both container implementations
                //particles.emplace_back(x_arg, v_arg, mass, generateNumber);
            }
        }
    }
   // container.setParticles(particles);
}


// TODO: How does the mesh width work for the sphere
void ParticleGenerator::instantiateSphere(ParticleContainer &container, std::array<double, 3> center, unsigned int nrMR, double h)
{
}
