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

public:
    SphereConstructor(std::array<double, 3> &cCoordinates, std::array<double, 3> &iVelocity,
                      int radius, double distance, double mass);

    const std::array<double, 3> getCenterCoordinates() const;

    const std::array<double, 3> getInitialVelocity() const;

    const int getRadius() const;

    const double getDistance() const;

    const double getMass() const;


};