#pragma once
#include "boundaryConditions/BoundaryCondition.h"

class Reflecting : public BoundaryCondition
{
public:
    Reflecting(double position, int direction, std::function<void(Particle &, Particle &)> &forceLambda, std::array<double, 3> domainSize);
	~Reflecting() = default;
	bool affectsForce() override;
	void applyBoundCondition(Particle &a) override;
    void applyHaloCondition(Particle &a) override;
};
