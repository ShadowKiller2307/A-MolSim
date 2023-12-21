
#include "CuboidConstructor.h"

CuboidConstructor::CuboidConstructor(std::array<double, 3>& llfc, std::array<unsigned int, 3>& particlesPerDimension,
                                     std::array<double, 3>& particleVelocity,
                                     double h, double mass, int type,double sigma, double epsilon) {
    this->llfc = llfc;
    this->particlesPerDimension = particlesPerDimension;
    this->particleVelocity = particleVelocity;
    this->h = h;
    this->mass = mass;
    this->type = type;
    this->sigma = sigma;
    this->epsilon =epsilon;

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

int CuboidConstructor::getType() const {
    return type;
}

double CuboidConstructor::getSigma() const {
    return sigma;
}
double CuboidConstructor::getEpsilon() const {
    return epsilon;
}