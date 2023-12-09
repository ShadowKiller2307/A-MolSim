#include "BoundaryCondition.h"
#include <iostream>

BoundaryCondition::BoundaryCondition(double position, int direction, std::function<void(Particle &, Particle &)> &forceLambda) : position_(position), direction_(direction), forceLambda_(forceLambda)
{
}

/*bool BoundaryCondition::affectsForce()
{
}*/

void BoundaryCondition::applyBoundCondition(Particle &a){}

bool BoundaryCondition::affectsForce() {
    std::cout << "Bitte mich nicht aufrufen!" <<std::endl;
    return false;
};
