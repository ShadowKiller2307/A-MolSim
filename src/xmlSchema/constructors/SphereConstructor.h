#pragma once

#include <iostream>
#include <array>

class SphereConstructor {

private:
    std::array<double, 3> centerCoordinates{};
    std::array<double, 3> initialVelocity{};
    int radius;
    double distance;
    double mass;
    double sigma;
    double epsilon;

public:
    SphereConstructor(std::array<double, 3> &cCoordinates, std::array<double, 3> &iVelocity,
                      int radius, double distance, double mass,double sigma,double epsilon);

    std::array<double, 3> getCenterCoordinates() const;

    std::array<double, 3> getInitialVelocity() const;

    int getRadius() const;

    double getDistance() const;

    double getMass() const;
    double getSigma() const;
    double getEpsilon() const;


};