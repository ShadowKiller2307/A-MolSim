#pragma once
#include <array>
#include <vector>
#include "ParticleContainerDS.h"
#include "Particle.h"
#include "utils/MaxwellBoltzmannDistribution.h"
#include "utils/ArrayUtils.h"
/**
 * @brief This class provide methods to create a cuboid or sphere of particles, defined by the command line or default values,
 * which will be added to a particleContainer
 */
class ParticleGenerator
{
public:
    /**
    * @brief instantiate a particle cuboid in the container
    * @param container container where the particles are added
    * @param llfc lower left front corner
    * @param particlePerDimension amount of particles in every dimension
    * @param particleVelocity initial velocity of the particles in the cuboid
    * @param h grid width
    * @param mass mass of the particles
    * @param generateNumber to mark to which Cuboid the particles belong (especially for containers which contain more than one cuboid)
    * @return void
    */
    void instantiateCuboid(ParticleContainer &container, std::array<double, 3> llfc, std::array<unsigned int, 3> particlePerDimension,
                           std::array<double, 3> particleVelocity, double h, double mass, int generateNumber);

    /**
     * instantiate a particle sphere in the container
     * @param container container where the particles are added
     * @param center the center of the sphere
     * @param nrMR number of molecules along the radius
     * @param h meshwidth between the molecules
     * @return void
     */
    void instantiateSphere(ParticleContainer &container, std::array<double, 3> center, unsigned int nrMR, double h);

    void instantiateCuboidNew(ParticleContainer &container, std::array<double, 3> llfc,
                                                 std::array<unsigned int, 3> particlePerDimension, std::array<double, 3> particleVelocity,
                                                 double h, double mass, int generateNumber);
};
