#pragma once
#include "boundaryConditions/BoundaryCondition.h"

class Outflow : public BoundaryCondition
{
public:
    Outflow(double position, int direction, std::function<void(Particle &, Particle &)> &forceLambda);
    ~Outflow() = default;
    bool affectsForce() override;
    void applyBoundCondition(Particle &a) override;
};
