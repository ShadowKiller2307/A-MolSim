#include "SphereConstructor.h"



SphereConstructor::SphereConstructor(std::array<double,3>& cCoordinates, std::array<double,3>& iVelocity,
int radius,double distance,double mass){
    this->centerCoordinates = cCoordinates;
    this->initialVelocity = iVelocity;
    this->radius = radius;
    this->distance = distance;
    this->mass = mass;
}


const std::array<double,3> SphereConstructor::getCenterCoordinates() const {
    return centerCoordinates;
}
const std::array<double,3> SphereConstructor::getInitialVelocity() const{
    return initialVelocity;
}

const int SphereConstructor::getRadius() const {
    return radius;
}

const double SphereConstructor::getDistance() const {
    return distance;
}

const double SphereConstructor::getMass() const {
    return mass;
}