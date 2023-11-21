/**
 * @brief This class implements the ParticleContainer using the linked cell algortihm
 */

#pragma once
#include "ParticleContainer.h"
#include "Particle.h"

using cell = std::vector<Particle>;

class ParticleContainerLC : ParticleContainer {
private:
    std::vector<cell> cells;
    std::array<double, 3> domainSize;
    double cutoffRadius;

public:
    void iterOverPairs(const std::function<void(Particle &a, Particle &b)> &f);
    /**
     * @brief iterate over the particles which are currectly located in the boundary zone
     */
    void iterBoundary();
    /**
     * @brief iterate over the particles which are currectly located in the halo zone
     */
     void iterHalo();

};
