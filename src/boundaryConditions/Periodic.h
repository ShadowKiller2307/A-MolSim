#pragma once
#include "BoundaryCondition.h"

class Periodic : public BoundaryCondition{
public:
    Periodic(double position, int direction, std::function<void(Particle &, Particle &)> &forceLambda, std::array<double, 3> domainSize);
    ~Periodic() = default;
    bool affectsForce() override;
    void applyBoundCondition(Particle &a, std::vector<cell> &particlesOtherSide) override;
    /**
     * @brief the particle should be moved across the boundary
     * @param a
     */
    void applyHaloCondition(Particle &a) override;
};