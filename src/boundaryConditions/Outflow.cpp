#include "boundaryConditions/Outflow.h"
#include "Outflow.h"

Outflow::Outflow(double position, int direction, std::function<void(Particle &, Particle &)> &forceLambda, std::array<double, 3> domainSize) : BoundaryCondition(position, direction, forceLambda, domainSize)
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

/**
 * @brief Particle a should be deleted
 * @param a
 */
void Outflow::applyHaloCondition(Particle &a) {
    BoundaryCondition::applyHaloCondition(a);
}

bool Outflow::affectsHalo() {
    return true;
}
