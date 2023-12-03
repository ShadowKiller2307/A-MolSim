#include "Reflecting.h"

Reflecting::Reflecting(double position, int direction, std::function<void(Particle &, Particle &)> &forceLambda) : BoundaryCondition(position, direction, forceLambda)
{
}

bool Reflecting::affectsForce()
{
	return true;
}

void Reflecting::applyBoundCondition(Particle &a)
{
	Particle ghostParticle = Particle();
	auto ghostPos = a.getX();
	ghostPos[direction_] = position_ + (position_ - ghostPos[direction_]);
	ghostParticle.setX(ghostPos);
	forceLambda_(a, ghostParticle);
}
