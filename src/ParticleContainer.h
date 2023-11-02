/**
 * @brief Class to encapsulate a vector of particles, it has an attribute ForceCalulator
 * to define which calculation method should be used on the particles
 */
#pragma once

#include "Particle.h"
#include <vector>

#include "ForceCalculator.h"

// Interne Implementierung die wie encapsulaten


class ParticleContainer {

private:
    std::vector<Particle> particles;
    ForceCalculator *forceCalculator;
    double deltaTTwo;


public:
    std::vector<Particle>* getParticles();

    /**
     * @brief calculate the force for all particles
     * @param None
     * @return void
    */
    void calculateForces();
    /**
     * @brief calculate the velocity for all particles
     * @param None
     * @return void
    */
    void calculateVelocity();
    /**
     * @brief calculate the position for all particles
     * @param None
     * @return void
    */
    void calculatePosition();


    explicit ParticleContainer();

    /**
     * @brief sets the particles for the container
     * @param particles1 the particle vector to be set for the container
     * @return void
    */
    void setParticles(const std::vector<Particle>& particles1);
    /**
     * @brief sets the ForceCalculator for the container
     * @param mode represent the different approaches for the force calculation
     * @return void
    */
    void setForceCalculator(int mode);

    /**
     * @brief sets the deltaT for the container
     * @param deltaT the deltaT passed by the user or the default value
     * @return void
    */
    void setDeltaTTwo(double deltaT);
};
