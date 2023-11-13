/**
 * @brief This class implements the approach of force calculation defined by Störmer-Verlet
 * it inherits the ForceCalculator class so that in can be passed to the ParticleContainer
 */

#pragma once

#include <vector>
#include "Particle.h"
#include "ForceCalculator.h"

class ForceV1 : public ForceCalculator
{

public:
    /**
     * @brief Implementation of the Störmer-Verlet approach for the forces
     * @param particles The particles from the ParticleContainer for which the force calculation will be executed
     * @return void
     */
    void calculateForces(std::vector<Particle> &particles) override;
    void calculateForcesWithLambda(std::vector<Particle> &particles) override;
};
