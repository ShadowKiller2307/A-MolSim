#pragma once

#include <vector>
#include "Particle.h"
#include <limits>

using cell = std::vector<Particle>;
class Thermostat{
private:
    const int NR_DIMENSIONS = 2;
    const double MAX_DOUBLE =std::numeric_limits<double>::max();
    double currentTemperature;
    double temperatureTarget;
    double maxTemperatureDifference;

    double initialTemperature;

public:
    //Thermostat(double initT,double tempTarget, double maxDiff = MAX_DOUBLE);
    Thermostat(double initT,double tempTarget, double maxDiff);
    void initializeTemperature(std::vector<cell> &cells);
    void scaleVelocity(std::vector<cell> &cells, double newTemp);
    void gradualScaleVelocity(std::vector<cell> &cells, double tempToAdd);
    void adjustTemperature();
    double calculateKinEnergy(std::vector<cell> &cells);
    double calculateCurrentTemperature(int nrParticles, double kineticEnergy);
    void regulateTemperature(std::vector<cell> &cells, int nrParticles);
};
