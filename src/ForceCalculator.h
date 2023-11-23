/**
 * @brief The ForceCalculator defines an abstract method calculateForces which is used for the different
 * force calculation approaches
 */
#pragma once

#include <vector>
#include "Particle.h"
#include "ParticleContainerDS.h"
#include "ParticleContainer.h"

class ParticleContainerDS; //needs to be defined so that the compiler doesn't throw an error

class ForceCalculator
{
public:
    /**
     * @brief virtual function which represents the different force calculation approaches
     * @param particles The particles from the ParticleContainerDS for which the force calculation will be executed
     * @return void
     */
    virtual void calculateForces(std::vector<Particle> &particles);
    /**
     * @brief virtual function which represents the different force calculation approaches with lambdas
     * @param particles The particles from the ParticleContainerDS for which the force calculation will be executed
     * @return void
     */
    virtual void calculateForcesWithLambda(ParticleContainer &container);
};