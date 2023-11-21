#include <vector>
#include "ParticleContainerDS.h"
#include "utils/ArrayUtils.h"
#include "HelperFunctions.h"
#include "ForceV1.h"
#include "LennardJonesForce.h"
#include <iostream>
#include "ParticleContainer.h"

// double deltaT{0.014};

ParticleContainerDS::ParticleContainerDS() = default;

std::vector<Particle> &ParticleContainerDS::getParticles() {
    return this->particles;
}

void ParticleContainerDS::setParticles(const std::vector<Particle> &particles1) {
    this->particles = particles1;
}

void ParticleContainerDS::calculateForces() {
    forceCalculator->calculateForces(particles);
}

void ParticleContainerDS::iterOverPairs(const std::function<void(Particle &a, Particle &b)> &forceLambda) {
    //debugLog("Currently applying iterOverPairs...\n");
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

    //debugLog("Currently applying calculatePosition...\n");
    int i = 0;
    for (auto &p: particles) {
        LogManager::getInstance().getLogger()->
                log(isDebug(), "Calculating position for particle number {}.\n", i);
        std::array<double, 3> force = p.getF();
        HelperFunctions::scalarOperations(force, 2 * p.getM(), true);
        HelperFunctions::scalarOperations(force, std::pow(deltaTTwo, 2), false);
        std::array<double, 3> newPosition = p.getX() + deltaTTwo * p.getV() + force;
        LogManager::getInstance().getLogger()
                ->log(isDebug(), "The new position for particle {} is {}.\n", i,
                      HelperFunctions::arrayToString(newPosition));
        p.setX(newPosition);

    }
}

void ParticleContainerDS::calculateVelocity() {
    LogManager::getInstance().getLogger()->log(isDebug(), "Currently applying calculateVelocity...\n");
    int i = 0;
    for (auto &p: particles) {
        LogManager::getInstance().getLogger()->
                log(isDebug(), "Calculating velocity for particle number {}.\n", i);
        double twoTimesMass = 2 * p.getM();
        std::array<double, 3> sumOfForces = p.getOldF() + p.getF();
        HelperFunctions::scalarOperations(sumOfForces, twoTimesMass, true);

        HelperFunctions::scalarOperations(sumOfForces, deltaTTwo, false);
        std::array<double, 3> newVelocity = p.getV() + sumOfForces;
        LogManager::getInstance().getLogger()
                ->log(isDebug(), "The new velocity for particle {} is {}.\n", i,
                      HelperFunctions::arrayToString(newVelocity));
        p.setV(newVelocity);
        i++;
    }
}

// https://sourcemaking.com/design_patterns/strategy/cpp/1 Looked here how the strategy pattern works
void ParticleContainerDS::setForceCalculator(int mode) {
    LogManager::getInstance().getLogger()->log(isDebug(), "Setting the mode for the forces to {}.", mode);
    if (mode == 0) {
        forceCalculator = new ForceV1();
    } else if (mode == 1) {
        forceCalculator = new LennardJonesForce{5, 1};
    } else {
        forceCalculator = new ForceV1(); // If the input mode isn't defined, the forceCalculator is set to ForceV1
    }
}

void ParticleContainerDS::setDeltaTTwo(double deltaT) {
    this->deltaTTwo = deltaT;
}

double ParticleContainerDS::getDeltaTwo() {
    return this->deltaTTwo;
}

spdlog::level::level_enum ParticleContainerDS::isDebug() {
    if (LogManager::getInstance().getLevel() == spdlog::level::debug) {
        return spdlog::level::debug;
    }
    return spdlog::level::off;
}
