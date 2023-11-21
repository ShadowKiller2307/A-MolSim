#pragma once
#include "ForceCalculator.h"
#include "ParticleContainerDS.h"

/**
 * @brief This class implements the approach of force calculation defined by the LennardJones force calculation
 * it inherits the ForceCalculator class so that in can be passed to the ParticleContainerDS
 */

class LennardJonesForce : public ForceCalculator
{
private:
    double epsilon = 5;
    double sigma = 1;

public:
    LennardJonesForce(double epsilon, double sigma);
    /**
     * @brief implementation of the LennardJones approach for the forces
     * @param particles The particles from the ParticleContainerDS for which the force calculation will be executed
     * @return void
     */
    void calculateForces(std::vector<Particle> &particles) override;
    /**
    * @brief implementation of the LennardJones approach for the forces with a lambda
    * @param container The container of the particles for which the force calculation will be executed
    * @return void
    */
    void calculateForcesWithLambda(ParticleContainerDS &container) override;
};
