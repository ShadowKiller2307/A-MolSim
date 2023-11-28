#include <functional>
#include "Particle.h"
#pragma once

class BoundaryCondition {

    /**
     * @brief it takes the boundary and the force lambda and returns a lambda
     * that takes on Particle and returns the according force
     */
public:
    virtual std::function<void(Particle &a)> applyBoundary(std::function<void(Particle &, Particle &)> &forceLambda, double position, int direction);

};

