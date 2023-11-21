#include "ParticleContainerLC.h"
#include "math.h"
using cell = std::vector<Particle>;
/*
ParticleContainerLC::ParticleContainerLC(std::array<double, 3> domainSize, double cutoffRadius) {
    this->domainSize = domainSize;
    this->cutoffRadius = cutoffRadius;
    // initialize the cells vector according to the domainSize and cutoff radius
    if (fmod(domainSize[0], cutoffRadius) == 0 && fmod(domainSize[1], cutoffRadius) == 0, fmod(domainSize[2], cutoffRadius) == 0) {
        double amountOfCells = (domainSize[0]/cutoffRadius) * (domainSize[1]/cutoffRadius) * (domainSize[2]/cutoffRadius);
        this->cells.reserve(amountOfCells);
        for(int i = 0; i < amountOfCells; ++i) { // intialize the cells
            cells.emplace_back(); //this should add a new cell
        }
    }
    else { // the domains size in each dimension has to be a multiple of the cutoff radius
        //TODO: maybe print error message
    }
}*/
