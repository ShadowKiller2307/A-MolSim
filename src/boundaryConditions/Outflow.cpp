#include "boundaryConditions/Outflow.h"
#include "Outflow.h"

Outflow::Outflow(double position, int direction, std::function<void(Particle &, Particle &)> &forceLambda) : BoundaryCondition(position, direction, forceLambda)
{
}

bool Outflow::affectsForce()
{
	return false;
}

void Outflow::applyBoundCondition(Particle &a)
{
	// do nothing, let the particle continue on its path
}
