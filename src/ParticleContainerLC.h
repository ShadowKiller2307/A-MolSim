/**
 * @brief This class implements the ParticleContainer using the linked cell algortihm
 */

#pragma once
#include "ParticleContainer.h"
#include "Particle.h"
#include "cmath"

class ParticleContainerLC : public ParticleContainer {

public:
    ParticleContainerLC(std::array<double, 3> domainSize, double cutoffRadius);

    /**
     * @brief iterate over the particles which are currectly located in the boundary zone
     */
     void iterBoundary();
    /**
     * @brief iterate over the particles which are currectly located in the halo zone
     */
     void iterHalo();

    void add(Particle &a) override;

    void iterOverPairs(const std::function<void(Particle &a, Particle &b)> &forceLambda) override;
};
