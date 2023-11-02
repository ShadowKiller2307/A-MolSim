/**
 * @brief The ForceCalculator defines an abstract method calculateForces which is used for the different
 * force calculation approaches
 */

#pragma once

class ForceCalculator{
public:
    virtual void calculateForces(std::vector<Particle>& particles) {};
};