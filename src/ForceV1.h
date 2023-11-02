/**
 * @brief This class implements the approach of force calculation defined by Störmer-Verlet
 * it inherits the ForceCalculator class so that in can be passed to the ParticleContainer
 */

#pragma once

#include <vector>
#include "Particle.h"
#include "ForceCalculator.h"

class ForceV1: public ForceCalculator{
public:
    void calculateForces(std::vector<Particle> &particles) override;
    //void calculateForces(std::vector<Particle> &particles);
};


