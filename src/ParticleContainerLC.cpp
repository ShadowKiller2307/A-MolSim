#include <variant>
#include "ParticleContainerLC.h"
#include "math.h"
using cell = std::vector<Particle>;

ParticleContainerLC::ParticleContainerLC(std::array<double, 3> domainSize, double cutoffRadius) {
    this->domainSize = domainSize;
    this->cutoffRadius = cutoffRadius;
    // initialize the cells vector according to the domainSize and cutoff radius
    if (fmod(domainSize[0], cutoffRadius) == 0 && fmod(domainSize[1], cutoffRadius) == 0, fmod(domainSize[2], cutoffRadius) == 0) {
        this->cells.reserve(amountOfCells);
        for(int i = 0; i < amountOfCells; ++i) { // intialize the cells
            cells.emplace_back(); //this should add a new cell
        }
        // TODO maybe initialise more cells for the halo area
        unsigned int amountOfHaloCells = 0;
        for(int i = 0; i < amountOfHaloCells; ++i) {
            cells.emplace_back();
        }
    }
    else { // the domains size in each dimension has to be a multiple of the cutoff radius
        //TODO: maybe print error message
    }
}

void ParticleContainerLC::add(Particle &a)  {
          // compute the cell to which the particle will be added
          double xIndex = trunc(a.getX()[0] / cutoffRadius);
          double yIndex = trunc(a.getX()[1]/cutoffRadius);
          double zIndex = trunc(a.getX()[2]/cutoffRadius);
          double index = xIndex + cellsX * yIndex + cellsX * cellsY * zIndex;
          cells.at(index).emplace_back(a);
}

//void ParticleContainer::iterOverPairs(const std::function<void(Particle &, Particle &)> &forceLambda) {}

void ParticleContainerLC::iterBoundary() {
    unsigned int acc = 0;
    for (auto cellTemp : cells) {
        acc += cellTemp.size();
    }
}

/**
 * @brief iterate over the pairs according to the linked cell algorithm
 * @param forceLambda
 */
void ParticleContainerLC::iterOverPairs(const std::function<void(Particle &, Particle &)> &forceLambda) {
    /*
     * iterate over the cells
     * if a new cell is iterated over, only calculate the forces in the cell itself
     * and the forces between the particles in the current cell and the particles on the right hand side of the cell and
     * over the current the cell
     */
    for (auto &p: particles) {
        auto oldForce = p.getF();
        std::array<double, 3> zero = {0.0, 0.0, 0.0};
        p.setF(zero);
        p.setOldF(oldForce);
    }
    for (int x = 0; x < amountOfCells; ++x) {
        cell current = cells.at(x);
        bool rightCell = x % cellsX == cellsX-1;
        bool upperCell = x < (amountOfCells-cellsX);
        bool rightUpperCell = rightCell and upperCell;
        /*
         * store the right and upper cells in a reference, if no right or upper cell exists
         * store the current cell in the variable(so that the variables always have a value)
         */
        cell &currentRight = rightCell ? cells.at(x+1) : cells.at(x);
        cell &currentUpper = upperCell ? cells.at(x + cellsX) : cells.at(x);
        cell &currentUpperRight = rightUpperCell ? cells.at(x + cellsX + 1) : cells.at(x);
        for (int i = 0; i < current.size(); ++i) {
            Particle &pi = current.at(i);
            // forces within the current cell
            for (int j = i + 1; j < current.size(); ++j) {
                Particle &pj = current.at(j);
                forceLambda(pi, pj);
            }
            if (rightCell) {
                for (auto & pj : currentRight) {
                    forceLambda(pi, pj);
                }
            }
            if (upperCell) {
                for (auto & pj : currentUpper) {
                    forceLambda(pi, pj);
                }
            }
            if (rightUpperCell) {
                for (auto & pj : currentUpperRight) {
                    forceLambda(pi, pj);
                }
            }
        }
    }
}



void ParticleContainerLC::iterHalo() {

}

/*void ParticleContainerLC::getSize() {

}*/
