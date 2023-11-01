
#pragma once

#include <vector>
#include "Particle.h"
#include "ForceCalculator.h"

class ForceV1: public ForceCalculator{
public:
    void calculateForces(std::vector<Particle> &particles) override;
    //void calculateForces(std::vector<Particle> &particles);
};


