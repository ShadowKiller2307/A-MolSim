#include "SphereConstructor.h"



SphereConstructor::SphereConstructor(std::array<double,3>& cCoordinates, std::array<double,3>& iVelocity,
int radius,double distance,double mass){
    this->centerCoordinates = cCoordinates;
    this->initialVelocity = iVelocity;
    this->radius = radius;
    this->distance = distance;
    this->mass = mass;
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