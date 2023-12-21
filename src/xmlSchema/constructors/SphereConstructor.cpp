#include "SphereConstructor.h"



SphereConstructor::SphereConstructor(std::array<double,3>& cCoordinates, std::array<double,3>& iVelocity,
int radius,double distance,double mass,double sigma,double epsilon,int type){
    this->centerCoordinates = cCoordinates;
    this->initialVelocity = iVelocity;
    this->radius = radius;
    this->distance = distance;
    this->mass = mass;
    this->sigma = sigma;
    this->epsilon = epsilon;
    this->type = type;
}


std::array<double,3> SphereConstructor::getCenterCoordinates() const {
    return centerCoordinates;
}
std::array<double,3> SphereConstructor::getInitialVelocity() const{
    return initialVelocity;
}

int SphereConstructor::getRadius() const {
    return radius;
}

double SphereConstructor::getDistance() const {
    return distance;
}

double SphereConstructor::getMass() const {
    return mass;
}
double SphereConstructor::getSigma() const {
    return sigma;
}
double SphereConstructor::getEpsilon() const {
    return epsilon;
}
int SphereConstructor::getType() const {
    return type;
}