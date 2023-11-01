#pragma once

#include "Particle.h"
#include <vector>

#include "ForceCalculator.h"

// Interne Implementierung die wie encapsulaten


class ParticleContainer {

private:
    std::vector<Particle> particles;
    static void scalarOperations(std::array<double,3> &array, double scalar, bool isDivision);
    double euclideanNorm(const std::array<double, 3> &arr);
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
