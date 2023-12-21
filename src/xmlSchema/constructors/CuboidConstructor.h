#pragma once
#include <iostream>
#include <array>
class CuboidConstructor{
private:
    std::array<double,3> llfc{};
    std::array<unsigned int,3> particlesPerDimension{};
    std::array<double,3> particleVelocity{};
    double h;
    double mass;
    int type;
    double sigma;
    double epsilon;

public:
    CuboidConstructor(std::array<double,3>& llfc,std::array<unsigned int,3>& particlesPerDimension,
                      std::array<double,3>& particleVelocity, double h, double mass,
                      int type,double sigma, double epsilon);

    // Getter functions

    const std::array<double, 3> &getLlfc() const;

    const std::array<unsigned int, 3> &getParticlesPerDimension() const;

    const std::array<double, 3> &getParticleVelocity() const;

    double getH() const;

    double getMass() const;

    int getType() const;
    double getSigma() const;
    double getEpsilon() const;
};
