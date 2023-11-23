
#include "CuboidConstructor.h"

CuboidConstructor::CuboidConstructor(std::array<double, 3>& llfc, std::array<unsigned int, 3>& particlesPerDimension,
                                     std::array<double, 3>& particleVelocity, double h, double mass, double type) {
    this->llfc = llfc;
    this->particlesPerDimension = particlesPerDimension;
    this->particleVelocity = particleVelocity;
    this->h = h;
    this->mass = mass;
    this->type = type;

}

const std::array<double, 3> &CuboidConstructor::getLlfc() const {
    return llfc;
}

const std::array<unsigned int, 3> &CuboidConstructor::getParticlesPerDimension() const {
    return particlesPerDimension;
}

const std::array<double, 3> &CuboidConstructor::getParticleVelocity() const {
    return particleVelocity;
}

double CuboidConstructor::getH() const {
    return h;
}

double CuboidConstructor::getMass() const {
    return mass;
}

double CuboidConstructor::getType() const {
    return type;
}
