#include "BoundaryCondition.h"

std::function<void(Particle &a)>
BoundaryCondition::applyBoundary(std::function<void(Particle &, Particle &)> &forceLambda, double position,
                                 int direction) {
    return std::function<void(Particle &)>();
}
