/**
 * @brief This class implements the approach of force calculation defined by Störmer-Verlet
 * it inherits the ForceCalculator class so that in can be passed to the ParticleContainerDS
 */

#pragma once

#include <vector>
#include "Particle.h"
#include "ForceCalculator.h"
#include "ParticleContainerDS.h"

class ForceV1 : public ForceCalculator
{

public:
    /**
     * @brief Implementation of the Störmer-Verlet approach for the forces
     * @param particles The particles from the ParticleContainerDS for which the force calculation will be executed
     * @return void
     */
    void calculateForces(std::vector<Particle> &particles) override;
    /**
    * @brief implementation of the Störmer-Verlet approach for the forces with a lambda
    * @param container The container of the particles for which the force calculation will be executed
    * @return void
    */
    void calculateForcesWithLambda(ParticleContainer &container) override;
};
