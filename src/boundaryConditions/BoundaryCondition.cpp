#include "BoundaryCondition.h"

BoundaryCondition::BoundaryCondition(double position, int direction, std::function<void(Particle &, Particle &)> &forceLambda) : position_(position), direction_(direction), forceLambda_(forceLambda)
{
}

bool BoundaryCondition::affectsForce()
{
	return false;
}

void BoundaryCondition::applyBoundCondition(Particle &a){};
