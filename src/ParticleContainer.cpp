#include <vector>
#include "ParticleContainer.h"
#include "utils/ArrayUtils.h"
#include "HelperFunctions.h"


double deltaT{0.014};

ParticleContainer::ParticleContainer() = default;

std::vector<Particle> ParticleContainer::getParticles() {
    return this->particles;
}

void ParticleContainer::setParticles(const std::vector<Particle>& particles1){
    this->particles = particles1;

}

void ParticleContainer::calculateForces(){
    forceCalculator.calculateForces(particles);
}

void ParticleContainer::calculatePosition() {

}

void ParticleContainer::calculateVelocity() {
    for (auto &p : particles) {
        double twoTimesMass = 2 * p.getM();
        std::array<double,3> sumOfForces = p.getOldF() + p.getF();
        HelperFunctions::scalarOperations(sumOfForces, twoTimesMass, true);

        HelperFunctions::scalarOperations(sumOfForces,deltaT,false);
        std::array<double, 3> newVelocity = p.getV() + sumOfForces;
        p.setV(newVelocity);
    }
}

/*
void ParticleContainer::scalarOperations(std::array<double,3> &array, double scalar, bool isDivision){
    if (isDivision) {
        for (double & i : array) {
            i /= scalar;
        }
    }
    else {
        for (double & i : array) {
            i *= scalar;
        }
    }
}

double ParticleContainer::euclideanNorm(const std::array<double, 3> &arr) {
    double sum = 0.0;
    for (const auto &element: arr) {
        sum += element * element;
    }
    return std::sqrt(sum);
}*/

void ParticleContainer::setForceCalculator(ForceCalculator &forceCalculator1) {
    this->forceCalculator = forceCalculator1;
}





