/**
 * @brief The ForceCalculator defines an abstract method calculateForces which is used for the different
 * force calculation approaches
 */
#pragma once

#include <vector>
#include "Particle.h"
#include "ParticleContainer.h"

class ParticleContainer; //needs to be defined so that the compiler doesn't throw an error

class ForceCalculator
{
public:
    /**
     * @brief virtual function which represents the different force calculation approaches
     * @param particles The particles from the ParticleContainer for which the force calculation will be executed
     * @return void
     */
    // virtual ~ForceCalculator() = default;
    void setParticleContainer(ParticleContainer &containerArg);
    virtual void calculateForces(std::vector<Particle> &particles);
    virtual void calculateForcesWithLambda(ParticleContainer &container);
};