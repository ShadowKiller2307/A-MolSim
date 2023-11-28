#pragma once
#include "boundaryConditions/BoundaryCondition.h"

class Reflecting : public BoundaryCondition{
    std::function<void(Particle &a)> applyBoundary(std::function<void(Particle &, Particle &)> &forceLambda, double position, int direction) override;
};

