#pragma once
#include <array>
#include <vector>
#include "ParticleContainer.h"
#include "Particle.h"
#include "utils/MaxwellBoltzmannDistribution.h"
#include "utils/ArrayUtils.h"

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
    void instantiateCuboid(ParticleContainer &container, std::array<double, 3> llfc, std::array<unsigned int, 3> particlePerDimension, double h, double mass,
                           std::array<double, 3> particleVelocity, int generateNumber);
};
