#include "particleContainers/ParticleContainerLinCel.h"
#include "boundaryConditions/BoundaryCondition.h"

ParticleContainerLinCel::ParticleContainerLinCel
(double deltaT, double endTime, std::array<double, 3> domainSize, double cutoffRadius, std::vector<BoundaryCondition> &conditions) : ParticleContainer(deltaT, endTime){
    domainSize_ = domainSize;
    cutoffRadius_ = cutoffRadius;
    conditions_ = conditions;
    // initialize the cells vector according to the domainSize and cutoff radius
    if (fmod(domainSize[0], cutoffRadius) == 0 && fmod(domainSize[1], cutoffRadius) == 0) {
        cells.reserve(amountOfCells);
        for(int i = 0; i < amountOfCells; ++i) { // intialize the cells
            cells.emplace_back(); //this should add a new cell
        }
    }
    else { // the domains size in each dimension has to be a multiple of the cutoff radius
        throw std::invalid_argument{"Domain dimensions have to be multiples of the cutoff radius"};
    }
}

void ParticleContainerLinCel::iterOverInnerPairs(const std::function<void(Particle & , Particle & )> &f) {
    /*
     * iterate over the cells
     * if a new cell is iterated over, only calculate the forces in the cell itself
     * and the forces between the particles in the current cell and the particles on the right hand side of the cell and
     * over the current the cell
     * Calc     > Calc
     *   ^     /
     *   |   /
     * Cell -- > Calc
     */
    for (int x = 0; x < amountOfCells; ++x) {
        if (x < cellsX || x >= (amountOfCells - cellsX) || x % cellsX == 0 || x % cellsX == (cellsX -1)) {
            continue; // skip the halo cells
        }
        cell current = cells.at(x);
        bool rightCell = (x % cellsX < cellsX-2);
        bool upperCell = x < (amountOfCells-2*cellsX);
        bool rightUpperCell = rightCell and upperCell;
        /*
         * store the right and upper cells in a reference, if no right or upper cell exists
         * store the current cell in the variable(so that the variables always have a value)
         */
        // TODO: maybe replace with pointers
        cell &currentRight = rightCell ? cells.at(x+1) : cells.at(x);
        cell &currentUpper = upperCell ? cells.at(x + cellsX) : cells.at(x);
        cell &currentUpperRight = rightUpperCell ? cells.at(x + cellsX + 1) : cells.at(x);
        for (int i = 0; i < current.size(); ++i) {
            Particle &pi = current.at(i);
            // forces within the current cell
            for (int j = i + 1; j < current.size(); ++j) {
                Particle &pj = current.at(j);
                f(pi, pj);
            }
            if (rightCell) {
                for (auto & pj : currentRight) {
                    f(pi, pj);
                }
            }
            if (upperCell) {
                for (auto & pj : currentUpper) {
                    f(pi, pj);
                }
            }
            if (rightUpperCell) {
                for (auto & pj : currentUpperRight) {
                    f(pi, pj);
                }
            }
        }
    }
}

void ParticleContainerLinCel::iterBoundary(std::array<const std::function<void(Particle &, Particle &)>, 4> &boundaryLambda) {
    int i = cellsX + 1;
    // from left to right
    for (i; i < cellsX-2; ++i) {  //lower row
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

void ParticleContainerLinCel::iterHalo() {
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

void ParticleContainerLinCel::add(const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg, double mass,
                                  int type) {
    // compute the cell to which the particle will be added
    double xIndex = trunc(x_arg[0] / cutoffRadius_) + 1 + cellsX; //The "+1 + cellsX" is added because of the halo cells
    double yIndex = trunc(x_arg[1]/cutoffRadius_) + 1 + cellsX;
    double index = xIndex + cellsX * yIndex;
    cells.at(index).emplace_back(x_arg, v_arg, mass, type);
}

ParticleContainerLinCel::~ParticleContainerLinCel() = default;
