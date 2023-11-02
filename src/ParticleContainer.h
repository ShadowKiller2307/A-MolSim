/**
 * @brief Class to encapsulate a vector of particles, it has an attribute ForceCalulator
 * to define which calculation method should used on the particles
 */
#pragma once

#include "Particle.h"
#include <vector>

#include "ForceCalculator.h"

// Interne Implementierung die wie encapsulaten


class ParticleContainer {

private:
    std::vector<Particle> particles;
    ForceCalculator forceCalculator;


public:
    std::vector<Particle> getParticles();

    void calculateForces();
    void calculateVelocity();
    void calculatePosition();


    explicit ParticleContainer();

    void setParticles(const std::vector<Particle>& particles1);
    void setForceCalculator(ForceCalculator& forceCalculator1);
};
