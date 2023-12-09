#pragma once
#include "boundaryConditions/BoundaryCondition.h"

class Reflecting : public BoundaryCondition
{
public:
	Reflecting(double position, int direction, std::function<void(Particle &, Particle &)> &forceLambda);
	~Reflecting() = default;
	bool affectsForce() override;
	void applyBoundCondition(Particle &a) override;
};
