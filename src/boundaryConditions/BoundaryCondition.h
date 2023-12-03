#include <functional>
#include "Particle.h"
#pragma once

class BoundaryCondition
{
protected:
    /// @brief the direction the normal of the boundary is facing, 0 = x, 1 = y, 2 = z
    int direction_;
    /// @brief position of the boundary along its axis
    double position_;
    /// @brief how the boundary should affect particles close to it
    std::function<void(Particle &, Particle &)> forceLambda_;

public:
    BoundaryCondition(double position, int direction, std::function<void(Particle &, Particle &)> &forceLambda);
    ~BoundaryCondition() = default;
    virtual void applyBoundCondition(Particle &a);
};
