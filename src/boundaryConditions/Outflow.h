#pragma once
#include "boundaryConditions/BoundaryCondition.h"

class Outflow : public BoundaryCondition
{
public:
    Outflow(double position, int direction, std::function<void(Particle &, Particle &)> &forceLambda, std::array<double, 3> domainSize);
    ~Outflow() = default;
    bool affectsForce() override;
    void applyBoundCondition(Particle &a) override;
    void applyHaloCondition(Particle &a) override;
};
