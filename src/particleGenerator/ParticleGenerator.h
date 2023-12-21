#pragma once
#include "forces/Force.h"
#include <stdint.h>
#include <vector>

class ParticleContainer;

/*
    to use, pass for example to a function
    optinals({.deltaT = 0.00069, .endTime = 7})
    to set h and mass but use te default values
    for meanVelocity and generateNumber
*/
struct SimParams
{
public:
    double deltaT = -1;
    double endTime = -1;
    std::string containerType = "";
    std::string boundaries = "";
    std::array<double, 3> domainSize = {-1, -1, -1};
    int writeFrequency = -1;
    double cutoffRadius = -1;
};

class particleGenerator
{
private:
    static int generateNumber_;
    static constexpr const double meanVelocity_ = 0.1;
    static constexpr const int h_ = 1.122462048;
    static constexpr const int m_ = 1.0;

public:
    /**
     * @brief instantiates a cuboid in the given ParticleContainer
     * @param container container in which the cuboid will be instantiated
     * @param llfc lower left front corner
     * @param particlesPerDimension amount of particles per dimension
     * @param particleVelocity initial velocity
     * @param h meshwidth
     * @param m mass of each particle
     * @param type
     * @return void
     */
    static void instantiateCuboid(ParticleContainer **container, const std::array<double, 3> &llfc,
                                  const std::array<unsigned int, 3> &particlesPerDimension, std::array<double, 3> &particleVelocity, double h,
                                  double m, int type, double epsilon, double sigma, double initT);

    /**
     * @brief instantiates a sphere in the given ParticleContainer
     * @param container container in which the cuboid will be instantiated
     * @param center center of the sphere
     * @param sphereRadius radius of the sphere
     * @param particleVelocity initial velocity
     * @param h meshwidth
     * @param m mass of each particle
     * @param is2D is it a 2D or 3D
     * @param type
     * @return void
     */
    static void instantiateSphere(ParticleContainer **container, const std::array<double, 3> &center,
                                  const int32_t &sphereRadius, const std::array<double, 3> &particleVelocity, double h, double m, bool is2D, int type, double epsilon, double sigma, double initT);

    static void instantiateJSON(ParticleContainer **container, const std::string &path, Force &force, SimParams clArgs);

    static void instantiatePicture(ParticleContainer **container, const std::string &path, Force &force, SimParams clArgs);

    static void instantiateTxt(ParticleContainer **container, const std::string &path, Force &force, SimParams clArgs);

    static void instantiateXML(ParticleContainer **container, std::string &path, Force &force, SimParams clArgs);
};
