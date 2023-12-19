#include "particleContainers/ParticleContainerLinCel.h"
#include "logOutputManager/LogManager.h"
#include <iostream>
#include <array>

/**
 * "rrrrrr"
 * "oooooo"
 * "pppppp"
 */

ParticleContainerLinCel::ParticleContainerLinCel(double deltaT, double endTime, int writeFrequency,
                                                 const std::array<double, 3> &domainSize,
                                                 const std::string &bounds, double cutoffRadius,
                                                 bool useThermostat, double nThermostat,
                                                 bool isGradual, double initT,
                                                 double tempTarget,
                                                 double maxDiff) : ParticleContainer(deltaT, endTime, writeFrequency)
{
    domainSize_ = domainSize;
    cutoffRadius_ = cutoffRadius;
    if (bounds.length() < 6)
    {
        LogManager::errorLog("Exactly six bounds have to be specified");
        exit(1);
    }
    for (int i = 0; i < 6; ++i)
    {
        const char c = bounds.at(i);
        // 0, 120, 0, 50, 1, 0 for WS3
        const double pos = i % 2 == 0 ? 0.0 : domainSize[i / 2];
        const int dir = i / 2;
        auto f = force_; // TODO: lennJon.innerPairs() at the moment, replace back with boundaryPairs later //.boundaryPairs();

        if (c == 'r')
        {
            conditions.push_back(BoundaryCondition::Reflecting);
        }
        else if (c == 'o')
        {
            conditions.push_back(BoundaryCondition::Outflow);
        }
        else if (c == 'p')
        {
            conditions.push_back(BoundaryCondition::Periodic);
        }
        else
        {
            LogManager::errorLog("Bound type {} is unknown", c);
            exit(1);
        }
    }
    /*
     * calculate the amount of cells the domain will be consisting of
     * ceil is used for the case the domainSize isn't a multiple of the cutoffRadius
     */
    if (domainSize_[0] < 0.0)
    {
        domainSize_[0] = -domainSize[0];
    }
    if (domainSize_[1] < 0.0)
    {
        domainSize_[1] = -domainSize_[1];
    }

    if (domainSize_[2] < 0.0)
    {
        domainSize_[2] = -domainSize_[2];
    }

    cellsX = static_cast<unsigned int>(ceil(domainSize_[0] / cutoffRadius_) + 2);
    cellsY = static_cast<unsigned int>(ceil(domainSize_[1] / cutoffRadius_) + 2);
    cellsZ = static_cast<unsigned int>(ceil(domainSize_[2] / cutoffRadius_) + 2);
    // amountOfCells = cellsX * cellsY + 2 * (cellsX + 2) + 2 * cellsY;
    amountOfCells = cellsX * cellsY * cellsZ;
    /**
     *  If the domain size isn't a multiple of the cutoff radius
        the last row and/or column of cells will only be a fraction of a cell
        for example domainSize{10, 9} with cutoff radius 3
        the last column of the cells will be only a (1 * cutoffRadius) cell instead of a (cutoffradius * cutoffRadius) Cell
     */
    for (uint32_t i = 0; i < amountOfCells; ++i)
    { // initialize the cells
        cells.emplace_back();
    }
    buildLookUp();
    // std::cout << "constructor end\n";
    // initialize the thermostat
    //thermostat{initT, tempTarget, maxDiff};
    //thermostat{initT, tempTarget, maxDiff};

}

void ParticleContainerLinCel::add(const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg, double mass, int type)
{
    if (x_arg[0] >= 0 && x_arg[0] <= domainSize_[0] && // in x bounds
        x_arg[1] >= 0 && x_arg[1] <= domainSize_[1] && // in y bounds
        x_arg[2] >= 0 && x_arg[2] <= domainSize_[2])   // in z bounds
    {
        // compute the cell to which the particle will be added
        cells.at(translate3DPosTo1D(x_arg)).emplace_back(x_arg, v_arg, mass, type);
    }
}

void ParticleContainerLinCel::simulateParticles()
{
    auto begin = std::chrono::high_resolution_clock::now();
    iteration_ = 0;

    while (startTime_ < endTime_)
    {
        if (iteration_ % 100 == 0)
        {
            std::vector<Particle> allParticles;
            allParticles.clear();
            for (auto &cell : cells)
            {
                for (auto &i : cell)
                {
                    allParticles.emplace_back(i);
                }
            }
            outManager_->plotParticles2(allParticles, iteration_);
            if (iteration_)
            {
                auto end = std::chrono::high_resolution_clock::now();
                size_t diff = std::chrono::duration_cast<std::chrono::seconds>(end - begin).count();
                auto remainig = (static_cast<double>(diff) / iteration_) * (endTime_ - startTime_) / deltaT_;
                std::string min, sec;
                if (remainig > 60)
                {
                    min = std::to_string(static_cast<int>(remainig / 60)) + "m, ";
                    remainig = std::fmod(remainig, 60);
                }
                else
                {
                    min = "0m, ";
                }
                sec = std::to_string(static_cast<int>(remainig)) + "s   \r";
                std::cout << "ETA: " << min << sec << std::flush;
            }
        }
        // check whether the Thermostat should be applied
       /* if (use(iteration_ % nThermostat) == 0) {
            //apply thermostat
        }*/
        // calculate new f
        calculateForces();
        // calculate new x
        calculatePosition();
        // calculate new v
        calculateVelocity();
        iteration_++;
        startTime_ += deltaT_;
    }
    //  TODO (ADD): Log
    auto end = std::chrono::high_resolution_clock::now();
    size_t diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    std::cout << "Output written, took " + std::to_string(diff) + " milliseconds. (about " + (iteration_ > diff ? std::to_string(iteration_ / diff) + " iter/ms" : std::to_string(diff / iteration_) + " ms/iter") + ") Terminating...\n";
    //TODO: calculate MUPS

}

void ParticleContainerLinCel::calculateForces()
{
    for (auto &cell : cells)
    {
        for (auto &p : cell)
        {
            auto oldForce = p.getF();
            std::array<double, 3> zero = {0.0, 0.0, 0.0};
            p.setF(zero);
            p.setOldF(oldForce);
        }
    }
    iterOverInnerPairs(force_);
    iterBoundary2();
}

void ParticleContainerLinCel::calculatePosition()
{
    std::vector<Particle> addBack;
    addBack.clear();
    for (uint32_t x = 0; x < cellsX; ++x)
    {
        for (uint32_t y = 0; y < cellsY; ++y)
        {
            for (uint32_t z = 0; z < cellsZ; ++z)
            {
                cell &currentCell = cells.at(translate3DIndTo1D(x, y, z));
                for (unsigned int i = 0; i < currentCell.size(); ++i)
                {
                    Particle &particle = currentCell.at(i);
                    std::array<double, 3> force = particle.getF();
                    double factor = std::pow(deltaT_, 2) / (2 * particle.getM());
                    force = factor * force;
                    std::array<double, 3> newPosition = (particle.getX()) + deltaT_ * (particle.getV()) + force; // calculate position
                    unsigned int afterCellIndex = translate3DPosTo1D(newPosition);
                    if (translate3DIndTo1D(x, y, z) == afterCellIndex)
                    {
                        particle.setX(newPosition);
                    }
                    else
                    {
                        // the particle has to move to a new cell
                        particle.setX(newPosition);
                        addBack.emplace_back(particle);
                        currentCell.erase(currentCell.begin() + i);
                        --i;
                    }
                }
            }
        }
    }
    // add all particles from addback to the cells
    for (auto &i : addBack)
    {
        if ((translate3DPosTo1D(i.getX())) >= cells.size() || (i.getX() != i.getX()))
        {
            continue;
        }
        cells.at(translate3DPosTo1D(i.getX())).emplace_back(i);
    }
    iterHalo();
}

void ParticleContainerLinCel::iterOverAllParticles(const std::function<void(ParticleContainerLinCel::cell::iterator)> &f)
{
    for (size_t x = 1; x < cellsX - 1; x++)
    {
        for (size_t y = 1; y < cellsY - 1; y++)
        {
            for (size_t z = 1; z < cellsZ - 1; z++)
            {
                cell &c = cells.at(translate3DIndTo1D(x, y, z));
                for (auto it = c.begin();;)
                {
                    if (it == c.end())
                    {
                        break;
                    }
                    f(it);
                    it++;
                }
            }
        }
    }
}

void ParticleContainerLinCel::iterOverInnerPairs(const std::function<void(Particle &, Particle &)> &f)
{
    // TODO: look at this for 3D
    for (uint32_t x = 1; x < cellsX - 1; ++x)
    {
        for (uint32_t y = 1; y < cellsY - 1; ++y)
        {
            cell &current = cells.at(translate3DIndTo1D(x, y, 1));
            for (unsigned int i = 0; i < current.size(); ++i)
            {
                // current cell
                Particle &ppi = current.at(i);
                for (size_t j = i + 1; j < current.size(); ++j)
                {
                    Particle &ppj = current.at(j);
                    calcF(ppi, ppj);
                }
                // right cell
                if (x != cellsX - 2)
                {
                    auto &rightCell = cells.at(translate3DIndTo1D(x + 1, y, 1));
                    for (auto &ppj : rightCell)
                    {
                        // check whether ppj is within the cutoff radius of ppi
                        if (ArrayUtils::L2Norm(ppj.getX() - ppi.getX()) <= cutoffRadius_)
                        {
                            calcF(ppi, ppj);
                        }
                    }
                }
                // upper cell
                if (y != cellsY - 2)

                {
                    auto &upperCell = cells.at(translate3DIndTo1D(x, y + 1, 1));
                    for (auto &ppj : upperCell)
                    {
                        // check whether ppj is within the cutoff radius of ppi
                        if (ArrayUtils::L2Norm(ppj.getX() - ppi.getX()) <= cutoffRadius_)
                        {
                            calcF(ppi, ppj);
                        }
                    }
                }
                // right upper cell
                if ((x != cellsX - 2) && (y != cellsY - 2))
                {
                    auto &rightUpperCell = cells.at(translate3DIndTo1D(x + 1, y + 1, 1));
                    for (auto &ppj : rightUpperCell)
                    {
                        // check whether ppj is within the cutoff radius of ppi
                        if (ArrayUtils::L2Norm(ppj.getX() - ppi.getX()) <= cutoffRadius_)
                        {
                            calcF(ppi, ppj);
                        }
                    }
                }
                // left upper cell
                if ((x > 1) && (y != cellsY - 2))
                {
                    auto &leftUpperCell = cells.at(translate3DIndTo1D(x - 1, y + 1, 1));
                    for (auto &ppj : leftUpperCell)
                    {
                        // check whether ppj is within the cutoff radius of ppi
                        if (ArrayUtils::L2Norm(ppj.getX() - ppi.getX()) <= cutoffRadius_)
                        {
                            calcF(ppi, ppj);
                        }
                    }
                }
            }
        }
    }
}

std::function<void(uint32_t x, uint32_t y, uint32_t z)> ParticleContainerLinCel::createPeriodicLambdaBoundary(int direction, int position)
{
    // TODO this should apply the forces from the opposite boundary cells to the particles in the current cell
    return [](uint32_t x, uint32_t y, uint32_t z) {
    };
}

std::function<void(uint32_t x, uint32_t y, uint32_t z)> ParticleContainerLinCel::createReflectingLambdaBoundary(int direction, int position)
{
    // TODO: change to reference
    auto lambda = [=](uint32_t x, uint32_t y, uint32_t z)
    {
        unsigned int cellIndex = translate3DIndTo1D(x, y, z);
        auto &cell = cells.at(cellIndex);
        for (auto &p : cell)
        {
            Particle ghostParticle = Particle();
            auto ghostPos = p.getX();
            ghostPos[direction] = position + (position - ghostPos[direction]);
            ghostParticle.setX(ghostPos);
            if (ArrayUtils::L2Norm(p.getX() - ghostParticle.getX()) <= std::pow(2, 1.0 / 6))
            {
                // TODO: another check whether the force is really repulsing
                calcF(p, ghostParticle);
            }
        }
    };
    return lambda;
}

void ParticleContainerLinCel::iterBoundary()
{
    // TODO: uncomment Z direction
    auto calculateBothPlanesInDirection = [&](uint32_t primaryDimension, uint32_t secondaryDimension1, uint32_t secondaryDimension2, int direction)
    {
        uint32_t i, j, k;
        uint32_t &x = direction == 0 ? i : j;
        uint32_t &y = direction == 1 ? i : direction == 0 ? j
                                                          : k;
        uint32_t &z = direction == 2 ? i : k;
        std::function<void(uint32_t x, uint32_t y, uint32_t z)> lambda;
        for (i = 1; i < primaryDimension - 1; i += (primaryDimension - 3))
        {
            switch (conditions[direction + (i == 1 ? 0 : 1)])
            {
            case BoundaryCondition::Outflow:
                lambda = [](uint32_t x, uint32_t y, uint32_t z) {};
                break;
            case BoundaryCondition::Reflecting:
                lambda = createReflectingLambdaBoundary(direction, domainSize_.at(direction + (i == 1 ? 0 : 1)));
                break;
            case BoundaryCondition::Periodic:
                break;
            }
            for (j = 1; j < secondaryDimension1 - 1; ++j)
            {
                for (k = 1; k < secondaryDimension2 - 1; ++k)
                {
                    lambda(x, y, z);
                }
            }
        }
    };
    calculateBothPlanesInDirection(cellsX, cellsY, cellsZ, 0);
    calculateBothPlanesInDirection(cellsY, cellsX, cellsZ, 1);
    // calculateBothPlanesInDirection(cellsZ, cellsX, cellsY, 2);
}

void ParticleContainerLinCel::iterBoundary2()
{
    uint32_t z = 1;
    if (conditions[2] == BoundaryCondition::Reflecting || conditions[2] == BoundaryCondition::Periodic)
    {
        for (uint32_t x = 1; x < cellsX - 1; ++x)
        {
            //std::cout << "Unten, Index: " << translate3DIndTo1D(x, 1, z) << std::endl;
            if (conditions[2] == BoundaryCondition::Reflecting) {
                // std::cout << "Passt!" << std::endl;
                auto lambda = createReflectingLambdaBoundary(1,0);
                lambda(x, 1, z);
                //  createRefectingForce(x, 1, z, 1, 0);
            }
            else
            {
                auto lambda = createPeriodicLambdaBoundary(1, 0);
                lambda(x, 1, z);
            }
        }
    }
    // from down right corner to top right corner
    if (conditions[1] == BoundaryCondition::Reflecting || conditions[1] == BoundaryCondition::Periodic)
    {
        for (uint32_t y = 1; y < cellsY - 1; ++y)
        {
            //std::cout << "Rechts, Index: " << translate3DIndTo1D(cellsX - 2, y, z) << std::endl;
            if (conditions[1] == BoundaryCondition::Reflecting) {
                auto lambda = createReflectingLambdaBoundary(0, static_cast<int>(domainSize_[0]));
                lambda(cellsX - 2, y, z);
                // createRefectingForce(cellsX - 2, y, z, 1, static_cast<int>(domainSize_[0]));
            }
            else
            {
                auto lambda = createPeriodicLambdaBoundary(1, 0);
                lambda(cellsX - 2, y, z);
            }
        }
    }
    if (conditions[3] == BoundaryCondition::Reflecting || conditions[3] == BoundaryCondition::Periodic)
    {
        // from top right corner to top left corner
        for (uint32_t x = cellsX - 2; x > 0; --x)
        {
            //std::cout <<"Oben, Index: " << translate3DIndTo1D(x, cellsY - 2, z) << std::endl;
            if (conditions[3] == BoundaryCondition::Reflecting) {
                auto lambda = createReflectingLambdaBoundary(1,static_cast<int> (domainSize_[1]));
                lambda(x, cellsY - 2, z);
                // createRefectingForce(x, cellsY - 2, z, 1, static_cast<int> (domainSize_[1]));
            }
            else
            {
                auto lambda = createPeriodicLambdaBoundary(1, static_cast<int>(domainSize_[1]));
                lambda(x, cellsY - 2, z);
            }
        }
    }
    // from top left corner to down left corner
    if (conditions[0] == BoundaryCondition::Reflecting || conditions[0] == BoundaryCondition::Periodic)
    {
        for (uint32_t y = cellsY - 2; y > 0; --y)
        {
            //    std::cout << "Links, Index: " << translate3DIndTo1D(1, y, z) << std::endl;
            if (conditions[0] == BoundaryCondition::Reflecting) {
                auto lambda = createReflectingLambdaBoundary(0,0);
                lambda(1, y, z);
                // createRefectingForce(1, y, z, 0,0);
            }
            else
            {
                auto lambda = createPeriodicLambdaBoundary(0, 0);
                lambda(1, y, z);
            }
        }
    }
}


std::function<void(uint32_t x, uint32_t y, uint32_t z)> ParticleContainerLinCel::createOutflowLambdaHalo()
{
    return [&](uint32_t x, uint32_t y, uint32_t z)
    {
        unsigned int _1DIndex = translate3DIndTo1D(x, y, z);
        cells.at(_1DIndex).clear();
    };
}

std::function<void(uint32_t x, uint32_t y, uint32_t z)> ParticleContainerLinCel::createPeriodicLambdaHalo(int direction, int position)
{
    // TODO
    return [](uint32_t x, uint32_t y, uint32_t z) {};
}

void ParticleContainerLinCel::iterHalo()
{
    // TODO: uncomment 3D
    auto calculateBothPlanesInDirection = [&](uint32_t primaryDimension, uint32_t secondaryDimension1, uint32_t secondaryDimension2, int direction)
    {
        uint32_t i, j, k;
        uint32_t &x = direction == 0 ? i : j;
        uint32_t &y = direction == 1 ? i : direction == 0 ? j
                                                          : k;
        uint32_t &z = direction == 2 ? i : k;
        std::function<void(uint32_t x, uint32_t y, uint32_t z)> lambda;
        for (i = 0; i < primaryDimension; i += (primaryDimension - 1))
        {
            switch (conditions[direction + (i == 1 ? 0 : 1)])
            {
            case BoundaryCondition::Outflow:
                lambda = createOutflowLambdaHalo();
                break;
            case BoundaryCondition::Reflecting:
                lambda = [](uint32_t x, uint32_t y, uint32_t z) {}; // do nothing
                break;
            case BoundaryCondition::Periodic:
                lambda /*move Lambda*/;
                break;
            }
            for (j = 0; j < secondaryDimension1; ++j)
            {
                for (k = 0; k < secondaryDimension2; ++k)
                {
                    lambda(x, y, z);
                }
            }
        }
    };
    calculateBothPlanesInDirection(cellsX, cellsY, cellsZ, 0);
    calculateBothPlanesInDirection(cellsY, cellsX, cellsZ, 1);
    // calculateBothPlanesInDirection(cellsZ, cellsX, cellsY, 2);
}

unsigned int ParticleContainerLinCel::translate3DIndTo1D(unsigned int x, unsigned int y, unsigned int z) const
{
    // return lookup.at(x).at(y).at(z);
    return (x + cellsX * y + cellsX * cellsY * z);
}

unsigned int ParticleContainerLinCel::translate3DPosTo1D(std::array<double, 3> position) const
{
    unsigned int xIndex;
    unsigned int yIndex;
    unsigned int zIndex;
    if (position[0] >= 0)
    {
        if (position[0] >= (domainSize_[0] + 1))
        {
            xIndex = cellsX - 1;
        }
        xIndex = static_cast<unsigned int>(floor(position[0] / cutoffRadius_)) + 1;
    }
    else
    {
        xIndex = 0;
    }
    if (position[1] >= 0)
    {
        if (position[1] >= (domainSize_[1] + 1))
        {
            yIndex = cellsY - 1;
        }
        yIndex = static_cast<unsigned int>(floor(position[1] / cutoffRadius_)) + 1;
    }
    else
    {
        yIndex = 0;
    }
    if (position[2] >= 0)
    {
        if (position[2] >= (domainSize_[2] + 1))
        {
            zIndex = cellsZ - 1;
        }
        zIndex = static_cast<unsigned int>(floor(position[2] / cutoffRadius_)) + 1;
    }
    else
    {
        zIndex = 0;
    }
    unsigned int index = xIndex + cellsX * yIndex + cellsX * cellsY * zIndex;
    return index;
}

void ParticleContainerLinCel::buildLookUp()
{
    lookup = std::vector<std::vector<std::vector<int>>>(cellsX, std::vector<std::vector<int>>(cellsY, std::vector<int>(cellsZ)));
    for (size_t x = 0; x < cellsX; x++)
    {
        for (size_t y = 0; y < cellsY; y++)
        {
            for (size_t z = 0; z < cellsZ; z++)
            {
                lookup[x][y][z] = x + cellsX * y + cellsX * cellsY * z;
            }
        }
    }
}

double ParticleContainerLinCel::calculateKinEnergy()
{
    double energy = 0.0;
    auto energyLambda = [&energy](ParticleContainerLinCel::cell::iterator it)
    {
        energy += (*it).getM() * std::inner_product((*it).getV().begin(), (*it).getV().end(), (*it).getV().begin(), 0.0) / 2;
    };
    iterOverAllParticles(energyLambda);
    return energy;
}

double ParticleContainerLinCel::calculateTemperature()
{
    // current simulations only in 2D
    auto numberofDimensions = 2;
    return (calculateKinEnergy() * 2) / (getAmountOfParticles() * numberofDimensions);
}

size_t ParticleContainerLinCel::getAmountOfParticles()
{
    size_t acc = 0;
    for (auto &cell : cells)
    {
        acc += cell.size();
    }
    return acc;
}

std::vector<Particle> ParticleContainerLinCel::getParticles()
{
    std::vector<Particle> returnVector;
    for (auto &cell : cells)
    {
        for (auto &i : cell)
        {
            returnVector.emplace_back(i);
        }
    }
    return returnVector;
}

unsigned int ParticleContainerLinCel::getAmountOfCells() const
{
    return amountOfCells;
}

std::vector<std::vector<Particle>> &ParticleContainerLinCel::getCells()
{
    return cells;
}
