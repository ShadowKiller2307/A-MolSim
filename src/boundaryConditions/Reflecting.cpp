//
// Created by alex on 28.11.23.
//

#include "Reflecting.h"

std::function<void(Particle &a)>
Reflecting::applyBoundary(std::function<void(Particle &, Particle &)> &forceLambda, double position, int direction) {
    return BoundaryCondition::applyBoundary(forceLambda, position, direction);
}
