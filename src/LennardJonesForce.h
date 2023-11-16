#pragma once
#include "ForceCalculator.h"
#include "ParticleContainer.h"

/**
 * @brief This class implements the approach of force calculation defined by the LennardJones force calculation
 * it inherits the ForceCalculator class so that in can be passed to the ParticleContainer
 */

class LennardJonesForce : public ForceCalculator
{
private:
    double epsilon = 5;
    double sigma = 1;

public:
    /**
     * @brief implementation of the Störmer-Verlet approach for the forces
     * @param particles The particles from the ParticleContainer for which the force calculation will be executed
     * @return void
     */
    LennardJonesForce(double epsilon, double sigma);
    void calculateForces(std::vector<Particle> &particles) override;
    void calculateForcesWithLambda(ParticleContainer &container) override;
};
