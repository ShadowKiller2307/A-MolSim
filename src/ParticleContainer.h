#pragma once

#include "Particle.h"
#include <vector>

#include "ForceCalculator.h"

// Interne Implementierung die wie encapsulaten


class ParticleContainer {

private:
    std::vector<Particle> particles;




public:
    std::vector<Particle> getParticles();
    void calculateForces(ForceCalculator calculator);
    void calculateVelocity();
    void calculatePosition();

    explicit ParticleContainer();

};
