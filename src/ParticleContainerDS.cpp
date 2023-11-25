#include <vector>
#include "ParticleContainerDS.h"
#include "utils/ArrayUtils.h"
#include "HelperFunctions.h"
#include "ForceV1.h"
#include "LennardJonesForce.h"
#include <iostream>
#include "ParticleContainer.h"


ParticleContainerDS::ParticleContainerDS() = default;

void ParticleContainerDS::add(Particle &a) {
    particles.emplace_back(a);
}

std::vector<Particle> &ParticleContainer::getParticles() {
    return this->particles;
}

void ParticleContainer::setParticles(const std::vector<Particle> &particles1) {
    this->particles = particles1;
}

void ParticleContainer::calculateForces() {
    forceCalculator->calculateForces(particles);
}

void ParticleContainerDS::iterOverPairs(const std::function<void(Particle &a, Particle &b)> &forceLambda) {

    ParticleContainer::debugLog("Currently applying iterOverPairs...\n");
    for (auto &p: particles) {
        auto oldForce = p.getF();
        std::array<double, 3> zero = {0.0, 0.0, 0.0};
        p.setF(zero);
        p.setOldF(oldForce);
    }
    for (size_t i = 0; i < particles.size() - 1; ++i) {
        Particle &pi = particles.at(i);
        for (size_t j = i + 1; j < particles.size(); ++j) {
            Particle &pj = particles.at(j);
            forceLambda(pi, pj);
        }
    }
}

void ParticleContainerDS::calculatePosition() {


    ParticleContainer::debugLog("Currently applying calculatePosition...\n");
    int i = 0;
    for (auto &p: particles) {

        ParticleContainer::debugLog("Calculating position for particle number {}.\n", i);
        std::array<double, 3> force = p.getF();
        HelperFunctions::scalarOperations(force, 2 * p.getM(), true);
        HelperFunctions::scalarOperations(force, std::pow(deltaTTwo, 2), false);
        std::array<double, 3> newPosition = p.getX() + deltaTTwo * p.getV() + force;
        ParticleContainer::debugLog("The new position for particle {} is {}.\n", i,
                                    HelperFunctions::arrayToString(newPosition));

        p.setX(newPosition);
        i++;
    }
}

void ParticleContainerDS::calculateVelocity() {

    ParticleContainer::debugLog("Currently applying calculateVelocity...\n");

    int i = 0;
    for (auto &p: particles) {

        ParticleContainer::debugLog("Calculating velocity for particle number {}.\n", i);
        double twoTimesMass = 2 * p.getM();
        std::array<double, 3> sumOfForces = p.getOldF() + p.getF();
        HelperFunctions::scalarOperations(sumOfForces, twoTimesMass, true);

        HelperFunctions::scalarOperations(sumOfForces, deltaTTwo, false);
        std::array<double, 3> newVelocity = p.getV() + sumOfForces;

        ParticleContainer::debugLog("The new velocity for particle {} is {}.\n", i,
                                    HelperFunctions::arrayToString(newVelocity));
        p.setV(newVelocity);
        i++;
    }
}



