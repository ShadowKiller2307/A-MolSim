#pragma once

#include <functional>
#include "Particle.h"
class BoundaryCondition
{
protected:

    /// @brief the direction the normal of the boundary is facing, 0 = x, 1 = y, 2 = z
    int direction_;
    /// @brief position of the boundary along its axis
    double position_;
    /// @brief how the boundary should affect particles close to it
    std::function<void(Particle &, Particle &)> forceLambda_;
    std::array<double, 3> domainSize_;

public:
    using cell = std::vector<Particle>;
    BoundaryCondition(double position, int direction, std::function<void(Particle &, Particle &)> &forceLambda, std::array<double, 3> domainSize);
    ~BoundaryCondition() = default;
    virtual bool affectsForce();
    virtual void applyHaloCondition(Particle &a);
    virtual void applyBoundCondition(Particle &a);
    virtual void applyBoundCondition(Particle &a, std::vector<cell> &particlesOtherSide);
};
