#include <variant>
#include "ParticleContainerLC.h"
#include "math.h"
#include "HelperFunctions.h"
#include "utils/ArrayUtils.h"
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

void ParticleContainerLC::calculateCellIndex(double xPos, double yPos, double zPos) {
    double xIndex = trunc(xPos / cutoffRadius);
    double yIndex = trunc(yPos/cutoffRadius);
    double zIndex = trunc(zPos/cutoffRadius);
    double index = xIndex + cellsX * yIndex + cellsX * cellsY * zIndex;
}

void ParticleContainerLC::add(Particle &a)  {
          // compute the cell to which the particle will be added
          double xIndex = trunc(a.getX()[0] / cutoffRadius);
          double yIndex = trunc(a.getX()[1]/cutoffRadius);
          double zIndex = trunc(a.getX()[2]/cutoffRadius);
          double index = xIndex + cellsX * yIndex + cellsX * cellsY * zIndex;
          cells.at(index).emplace_back(a);
}

void ParticleContainerLC::iterBoundary(std::array <const std::function<void(Particle &, Particle &)>, 4> &boundaryLambda) {
    int i = cellsX + 1;
    // from left to right
    for (i; i < cellsX-2; ++i) {
        cell &currentCell = cells.at(i);
    }
    // from down to up
    for (int j = 0; j < cellsY - 3; ++j) {
        i += cellsX;
        cell &currentCell = cells.at(i);
    }
    //from right to left
    for (int j = 0; j < cellsX-3; ++j) {
        --i;
        cell &currentCell = cells.at(i);
    }
    // from up to down
    for (int j = 0; j < cellsY-4; ++j) {
        i = i - cellsX;
        cell &currentCell = cells.at(i);
    }
}


// iterate over the halo cells like this
//         <-<-<-^
//         |     |
//  Start ->->->->
/**
 * @brief Clear the particles of every halo cell
 */
void ParticleContainerLC::iterHalo() {
    int i = 0;
    // from left to right
    for (i; i < cellsX; ++i) {
        cell &currentCell = cells.at(i);
        currentCell.clear(); // delete all particles from the current halo cell
    }
    // from down to up
    for (int j = 0; j < cellsY - 1; ++j) {
        i += cellsX;
        cell &currentCell = cells.at(i);
        currentCell.clear();
    }
    //from right to left
    for (int j = 0; j < cellsX-1; ++j) {
        --i;
        cell &currentCell = cells.at(i);
        currentCell.clear();
    }
    // from up to down
    for (int j = 0; j < cellsY-2; ++j) {
        i = i - cellsX;
        cell &currentCell = cells.at(i);
        currentCell.clear();
    }
}

/**
 * @brief iterate over the pairs according to the linked cell algorithm
 * @param forceLambda
 */
void ParticleContainerLC::iterOverPairs(const std::function<void(Particle &, Particle &)> &forceLambda) {
    //TODO: maybe dont iterate over the halo cells here?
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

void ParticleContainerLC::calculatePosition() {
    LogManager::debugLog("Currently applying calculatePosition...\n");
    int i = 0;
    for(int j = 0; j < cells.size(); ++j) {
    //for (auto &currentCell: cells) {
        auto &currentCell = cells.at(j);
        for(int x = 0; x < currentCell.size(); ++x) {
        //for (auto &p : currentCell) {
            auto &p = currentCell.at(x);
            LogManager::debugLog("Calculating position for particle number {}.\n", i);
            std::array<double, 3> force = p.getF();
            HelperFunctions::scalarOperations(force, 2 * p.getM(), true);
            HelperFunctions::scalarOperations(force, std::pow(deltaTTwo, 2), false);
            std::array<double, 3> newPosition = p.getX() + deltaTTwo * p.getV() + force;
            LogManager::debugLog("The new position for particle {} is {}.\n", i,
                                        HelperFunctions::arrayToString(newPosition));
            //update the cell
            /*if (newPosition[0] < 0 || newPosition[0] > domainSize[0] || newPosition[1] < 0 || newPosition[1] > domainSize[1]) {
                //add the particle to the halo cell or delete it
            }*/
            //check whether the particle left the current cell
            double index = trunc(newPosition[0] / cutoffRadius) //TODO: could maybe be implemented in a more performant way
                    + cellsX * trunc(newPosition[1]/cutoffRadius); +
                    + cellsX * cellsY * trunc(newPosition[2]/cutoffRadius);;
            if (j != index) { // Particle has to be deleted in the old cell and added to the new cell
                currentCell.erase(currentCell.begin() + x);
                cells.at(index).emplace_back(p);
            }
            p.setX(newPosition);
            i++;
        }
    }
}

void ParticleContainerLC::calculateVelocity() {
    LogManager::debugLog("Currently applying calculateVelocity...\n");
    int i = 0;
    for (auto &currentCell: cells) {
        for (auto &p : currentCell) {

            LogManager::debugLog("Calculating velocity for particle number {}.\n", i);
            double twoTimesMass = 2 * p.getM();
            std::array<double, 3> sumOfForces = p.getOldF() + p.getF();
            HelperFunctions::scalarOperations(sumOfForces, twoTimesMass, true);

            HelperFunctions::scalarOperations(sumOfForces, deltaTTwo, false);
            std::array<double, 3> newVelocity = p.getV() + sumOfForces;

            LogManager::debugLog("The new velocity for particle {} is {}.\n", i,
                                        HelperFunctions::arrayToString(newVelocity));
            p.setV(newVelocity);
            i++;
        }
    }
}

/*void ParticleContainerLC::getSize() {

}*/
