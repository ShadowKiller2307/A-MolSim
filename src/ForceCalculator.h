/**
 * @brief The ForceCalculator defines an abstract method calculateForces which is used for the different
 * force calculation approaches
 */

#pragma once

class ForceCalculator{
public:
    /**
     * @brief virtual function which represents the different force calculation approaches
     * @param particles The particles from the ParticleContainer for which the force calculation will be executed
     * @return void
     */
    virtual void calculateForces(std::vector<Particle>& particles) {};
};