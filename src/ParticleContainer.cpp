#include <vector>
#include "ParticleContainer.h"
#include "utils/ArrayUtils.h"
#include "HelperFunctions.h"
#include "ForceV1.h"
#include "LennardJonesForce.h"
#include <iostream>

// double deltaT{0.014};

ParticleContainer::ParticleContainer() = default;

std::vector<Particle> &ParticleContainer::getParticles()
{
    return this->particles;
}

void ParticleContainer::setParticles(const std::vector<Particle> &particles1)
{
    this->particles = particles1;
}

void ParticleContainer::calculateForces()
{
    forceCalculator->calculateForces(particles);
}

void ParticleContainer::iterOverPairs(const std::function<void(Particle &a, Particle &b)> &forceLambda)
{
    for (auto &p : particles)
    {
        auto oldForce = p.getF();
        std::array<double, 3> zero = {0.0, 0.0, 0.0};
        p.setF(zero);
        p.setOldF(oldForce);
    }
    for (size_t i = 0; i < particles.size() - 1; ++i)
    {
        Particle &pi = particles.at(i);
        for (size_t j = i + 1; j < particles.size(); ++j)
        {
            Particle &pj = particles.at(j);
            forceLambda(pi, pj);
        }
    }
}

void ParticleContainer::calculatePosition()
{
    for (auto &p : particles)
    {
        std::array<double, 3> force = p.getF();
        HelperFunctions::scalarOperations(force, 2 * p.getM(), true);
        HelperFunctions::scalarOperations(force, std::pow(deltaTTwo, 2), false);
        std::array<double, 3> newPosition = p.getX() + deltaTTwo * p.getV() + force;
        p.setX(newPosition);
    }
}

void ParticleContainer::calculateVelocity()
{
    for (auto &p : particles)
    {
        double twoTimesMass = 2 * p.getM();
        std::array<double, 3> sumOfForces = p.getOldF() + p.getF();
        HelperFunctions::scalarOperations(sumOfForces, twoTimesMass, true);

        HelperFunctions::scalarOperations(sumOfForces, deltaTTwo, false);
        std::array<double, 3> newVelocity = p.getV() + sumOfForces;
        p.setV(newVelocity);
    }
}

// https://sourcemaking.com/design_patterns/strategy/cpp/1 Looked here how the strategy pattern works
void ParticleContainer::setForceCalculator(int mode)
{
    if (mode == 0)
    {
        forceCalculator = new ForceV1();
    }
    else if (mode == 1)
    {
        forceCalculator = new LennardJonesForce{5, 1};
    }
    else
    {
        forceCalculator = new ForceV1(); // If the input mode isn't defined, the forceCalculator is set to ForceV1
    }
}

void ParticleContainer::setDeltaTTwo(double deltaT)
{
    this->deltaTTwo = deltaT;
}

double ParticleContainer::getDeltaTwo()
{
    return this->deltaTTwo;
}
