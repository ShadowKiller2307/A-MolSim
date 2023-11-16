#pragma once
#include <array>
#include <vector>
#include "ParticleContainer.h"
#include "Particle.h"
#include "utils/MaxwellBoltzmannDistribution.h"
#include "utils/ArrayUtils.h"
/**
 * @brief This class creates a cuboid of particles, definded by the command line or default values,
 * which will be added to a particleContainer
 */
class ParticleGenerator
{
private:
    /*std::array<double, 3> llfc; // lower left frontside corner
    std::array<double, 3> particlePerDimension;
    double h;
    double mass;
    std::array<double, 3> particleVelocity;*/

    // TODO: last parameter can be hardcoded but maybe add later

    // TODO: Hier werden wir einen Cuboid initialisieren
public:
    void instantiateCuboid(ParticleContainer &container, std::array<double, 3> llfc, std::array<unsigned int, 3> particlePerDimension,
                           std::array<double, 3> particleVelocity, double h, double mass, int generateNumber);
};
