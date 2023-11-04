#pragma once
#include <array>

class Cuboid {
private:
    std::array<double, 3> llfc; // lower left frontside corner
    std::array<double, 3> particlePerDimension;
    double h;
    double mass;
    std::array<double, 3> particleVelocity;
    //TODO: last parameter can be hardcoded but maybe add later

};