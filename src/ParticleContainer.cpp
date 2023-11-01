#include <iostream>
#include <vector>
#include "ParticleContainer.h"

ParticleContainer::ParticleContainer() {
}

std::vector<Particle> ParticleContainer::getParticles() {
    return this->particles;
}

void ParticleContainer::calculateForces(ForceCalculator calculator){
    calculator.calculateForces();
}

void ParticleContainer::calculatePosition() {

}

void ParticleContainer::calculateVelocity() {

}