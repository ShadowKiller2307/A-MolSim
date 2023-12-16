#include "BoundaryCondition.h"
#include <iostream>

BoundaryCondition::BoundaryCondition(double position, int direction, std::function<void(Particle &, Particle &)> &forceLambda, std::array<double, 3> domainSize) : position_(position), direction_(direction), forceLambda_(forceLambda),
domainSize_(domainSize)
{
}

/*bool BoundaryCondition::affectsForce()
{
}*/

void BoundaryCondition::applyBoundCondition(Particle &a){}

bool BoundaryCondition::affectsForce() {
    std::cout << "Bitte mich nicht aufrufen!" <<std::endl;
    return false;
}

void BoundaryCondition::applyHaloCondition(Particle &a) {
}

void BoundaryCondition::applyBoundCondition(Particle &a, std::vector<cell> &particlesOtherSide) {

}

