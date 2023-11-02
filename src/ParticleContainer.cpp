#include <vector>
#include "ParticleContainer.h"
#include "utils/ArrayUtils.h"
#include "HelperFunctions.h"


//double deltaT{0.014};

ParticleContainer::ParticleContainer() = default;

std::vector<Particle>* ParticleContainer::getParticles() {
    return &this->particles;
}

void ParticleContainer::setParticles(const std::vector<Particle> &particles1) {
    this->particles = particles1;

}

void ParticleContainer::calculateForces() {
    forceCalculator.calculateForces(particles);
}

void ParticleContainer::calculatePosition() {
    for (auto &p: particles) {
        std::array<double, 3> force = p.getF();
        HelperFunctions::scalarOperations(force, 2 * p.getM(), true);
        HelperFunctions::scalarOperations(force, std::pow(deltaTTwo, 2), false);
        std::array<double, 3> newPosition = p.getX() + deltaTTwo * p.getV() + force;
        p.setX(newPosition);
    }
}

void ParticleContainer::calculateVelocity() {
    for (auto &p: particles) {
        double twoTimesMass = 2 * p.getM();
        std::array<double, 3> sumOfForces = p.getOldF() + p.getF();
        HelperFunctions::scalarOperations(sumOfForces, twoTimesMass, true);

        HelperFunctions::scalarOperations(sumOfForces, deltaTTwo, false);
        std::array<double, 3> newVelocity = p.getV() + sumOfForces;
        p.setV(newVelocity);
    }
}

void ParticleContainer::setForceCalculator(ForceCalculator &forceCalculator1) {
    this->forceCalculator = forceCalculator1;
}

void ParticleContainer::setDeltaTTwo(double deltaT) {
    this->deltaTTwo = deltaT;
}





