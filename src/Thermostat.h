#pragma once

#include <vector>
#include "Particle.h"
#include <limits>

using cell = std::vector<Particle>;
class Thermostat{
private:
    static const int NR_DIMENSIONS = 2;
    constexpr static const double MAX_DOUBLE =std::numeric_limits<double>::max();
    static double currentTemperature;
    double temperatureTarget;
    static double maxTemperatureDifference;

    double initialTemperature;

public:
    Thermostat(double initT,double tempTarget, double maxDiff = MAX_DOUBLE);
    void initializeTemperature(std::vector<cell> &cells);
    static void scaleVelocity(std::vector<cell> &cells, double newTemp);
    static void gradualScaleVelocity(std::vector<cell> &cells, double tempToAdd);
    void adjustTemperature();
    static double calculateKinEnergy(std::vector<cell> &cells);
    static  double calculateCurrentTemperature(int nrParticles, double kineticEnergy);
    void regulateTemperature(std::vector<cell> &cells, int nrParticles);
};
