#pragma once
#include "boundaryConditions/BoundaryCondition.h"

class Overflow : public BoundaryCondition{
    std::function<void(Particle &a)> applyBoundary(std::function<void(Particle &, Particle &)> &forceLambda, double position, int direction) override;

};

