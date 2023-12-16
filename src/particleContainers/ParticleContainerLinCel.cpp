#include "particleContainers/ParticleContainerLinCel.h"
#include "boundaryConditions/Reflecting.h"
#include "boundaryConditions/Outflow.h"
#include "logOutputManager/LogManager.h"
#include <iostream>
#include <array>
#include <memory>

/**
 * "rrrrrr"
 * "oooooo"
 * "pppppp"
 */

/*enum class BoundaryCondition2 {
    Reflecting2,
    Periodic2,
    Outflow2
};*/

ParticleContainerLinCel::ParticleContainerLinCel(double deltaT, double endTime, int writeFrequency, const std::array<double, 3> &domainSize, const std::string &bounds, Force &force, double cutoffRadius) : ParticleContainer(deltaT, endTime, writeFrequency, force.innerPairs())
{
    //std::cout << "Constructor begin\n";
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
        auto f = force.boundaryPairs();
        // std::unique_ptr<BoundaryCondition> temp;
        if (c == 'r')
        {
           // reflectingBounds.emplace_back(pos, dir, f);
            conditions_.emplace_back(std::make_unique<Reflecting>((reflectingBounds.emplace_back(pos, dir, f, domainSize_))));
        }
        else if (c == 'o')
        {
            conditions_.emplace_back(std::make_unique<Outflow>((outflowBounds.emplace_back(pos, dir, f, domainSize_))));
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

    if (fmod(domainSize_[0], cutoffRadius_) != 0.0)
    { // last column of cells are smaller in width
        rightModulo = true;
    }
    if (fmod(domainSize_[1], cutoffRadius_) != 0.0)
    { // last row of cells are smaller in height
        upperModulo = true;
    }
    if (fmod(domainSize_[2], cutoffRadius_) != 0.0)
    {
        depthModulo = true;
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
    cellPointerNeedUpdate = true;
    //std::cout << "constructor end\n";
}

void ParticleContainerLinCel::iterOverInnerPairs(const std::function<void(Particle &, Particle &)> &f)
{
    //std::cout << "iterOverInnerPairs begin\n";
    for (uint32_t x = 1; x < cellsX - 1; ++x)
    {
        for (uint32_t y = 1; y < cellsY - 1; ++y)
        {
            for (uint32_t z = 1; z < cellsZ - 1; ++z)
            {
                cell &current = cells.at(translate3DIndTo1D(x, y, z));
                for (auto ppi : current)
                {
                    // right cell
                    if (x != cellsX - 2)
                    {
                        for (auto ppj : cells.at(translate3DIndTo1D(x + 1, y, z)))
                        {
                            // check whether ppj is within the cutoff radius of ppi
                            if (ArrayUtils::L2Norm(ppj.getX() - ppi.getX()) <= cutoffRadius_)
                            {
                                f(ppi, ppj);
                            }
                        }
                    }
                    // upper cell
                    if (y != cellsY - 2)
                    {
                        for (auto ppj : cells.at(translate3DIndTo1D(x, y + 1, z)))
                        {
                            // check whether ppj is within the cutoff radius of ppi
                            if (ArrayUtils::L2Norm(ppj.getX() - ppi.getX()) <= cutoffRadius_)
                            {
                                f(ppi, ppj);
                            }
                        }
                    }
                    // right upper cell
                    if ((x != cellsX - 2) && (y != cellsY - 2))
                    {
                        for (auto ppj : cells.at(translate3DIndTo1D(x + 1, y + 1, z)))
                        {
                            // check whether ppj is within the cutoff radius of ppi
                            if (ArrayUtils::L2Norm(ppj.getX() - ppi.getX()) <= cutoffRadius_)
                            {
                                f(ppi, ppj);
                            }
                        }
                    }
                    // left upper cell
                    if ((x > 1) && (y != cellsY - 2))
                    {
                        for (auto ppj : cells.at(translate3DIndTo1D(x - 1, y + 1, z)))
                        {
                            // check whether ppj is within the cutoff radius of ppi
                            if (ArrayUtils::L2Norm(ppj.getX() - ppi.getX()) <= cutoffRadius_)
                            {
                                f(ppi, ppj);
                            }
                        }
                    }
                    // TODO: expand implementation for 3D cases

                    // z + 1

                    // right cell

                    // upper cell

                    // right upper cell

                    // left upper cell
                }
            }
        }
    }
    //std::cout << "iterOverInnerPairs end\n";
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
                }
            }
        }
    }
}

/**
 * order:
 * left plane
 * right plane
 * lower plane
 * upper plane
 * back plane
 * front plane
 */
void ParticleContainerLinCel::iterBoundary()
{
    //std::cout << "iterBoundaries begin\n";
    for (uint32_t z = 1; z < cellsZ - 1; ++z)
    {
        if (z == 1)
        {
            if (conditions_[4]->affectsForce())
            {
                for (uint32_t x = 1; x < cellsX - 1; ++x)
                {
                    for (uint32_t y = 1; y < cellsY - 1; ++y)
                    {
                        cell &currentCell = cells.at(translate3DIndTo1D(x, y, z));
                        for (auto ppi : currentCell)
                        {
                            conditions_[4]->applyBoundCondition(ppi);
                        }
                    }
                }
            }
        }
        if (z == cellsZ - 2)
        {
            if (conditions_[5]->affectsForce()) {
                for (uint32_t x = 1; x < cellsX - 1; ++x) {
                    for (uint32_t y = 1; y < cellsY - 1; ++y) {
                        cell &currentCell = cells.at(translate3DIndTo1D(x, y, z));
                        for (auto ppi: currentCell) {
                            conditions_[5]->applyBoundCondition(ppi);
                        }
                    }
                }
            }
        }
        // iterate over the outer ring of inner cells for each inner z index
        //  from down left corner to down right corner
        if (conditions_[2]->affectsForce())
        {
            for (uint32_t x = 1; x < cellsX - 1; ++x)
            {
                cell &currentCell = cells.at(translate3DIndTo1D(x, 1, z));
                for (auto ppi : currentCell)
                {
                    conditions_[2]->applyBoundCondition(ppi);
                }
            }
        }
        // from down right corner to top right corner
        if (conditions_[1]->affectsForce())
        {
            for (uint32_t y = 1; y < cellsY - 1; ++y)
            {
                cell &currentCell = cells.at(translate3DIndTo1D(cellsX - 2, y, z));
                for (auto ppi : currentCell)
                {
                    conditions_[1]->applyBoundCondition(ppi);
                }
            }
        }
        if (conditions_[3]->affectsForce())
        {
            // from top right corner to top left corner

            for (uint32_t x = cellsX - 2; x > 0; --x)
            {
                cell &currentCell = cells.at(translate3DIndTo1D(x, cellsY - 2, z));
                for (auto ppi : currentCell)
                {
                    conditions_[3]->applyBoundCondition(ppi);
                }
            }
        }
        // from top left corner to down left corner
        if (conditions_[0]->affectsForce())
        {
            for (uint32_t y = cellsY - 2; y > 0; --y)
            {
                cell &currentCell = cells.at(translate3DIndTo1D(1, y, z));
                for (auto &ppi : currentCell)
                {
                    conditions_[0]->applyBoundCondition(ppi);
                }
            }
        }
    }
    //std::cout << "iterBoundaries end\n";
}

void ParticleContainerLinCel::iterHalo() {
    for (uint32_t x = 0; x < cellsX; ++x)
    {
        for (uint32_t y = 0; y < cellsY; ++y)
        {
            for (uint32_t z = 0; z < cellsZ; ++z)
            {
                if ((x == 0) || (x == cellsX - 1) || (y == 0) || (y == cellsY - 1) || (z == 0) || (z == cellsZ - 1))
                {
                    cell &currentCell = cells.at(translate3DIndTo1D(x, y, z));
                    currentCell.clear();
                }
            }
        }
    }

}



void ParticleContainerLinCel::add(const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg, double mass, int type)
{
    if (x_arg[0] >= 0 && x_arg[0] <= domainSize_[0] && // in x bounds
        x_arg[1] >= 0 && x_arg[1] <= domainSize_[1] && // in y bounds
        x_arg[2] >= 0 && x_arg[2] <= domainSize_[2])   // in z bounds
    {
        /*auto particle = new Particle(x_arg, v_arg, mass, type);
        // add the particle
        particles_.push_back(particle);*/
        // compute the cell to which the particle will be added
        cells.at(translate3DPosTo1D(x_arg)).emplace_back(x_arg, v_arg, mass, type);
    }
}

void ParticleContainerLinCel::calculatePosition() {
   // std::cout << "Richtiges calculatePosition\n";
    std::vector<Particle> addBack;
    for (uint32_t x = 0; x < cellsX; ++x)
    {
        for (uint32_t y = 0; y < cellsY; ++y)
        {
            for (uint32_t z = 0; z < cellsZ; ++z)
            {
               cell &currentCell = cells.at(translate3DIndTo1D(x, y, z));
               std::vector<unsigned int> indexesToDelete;
                for (unsigned int i = 0; i < currentCell.size(); ++i) {
                    Particle &particle = currentCell.at(i);
                    std::array<double, 3> force = particle.getF();
                    double factor = std::pow(deltaT_, 2) / (2 * particle.getM());
                    force = factor * force;
                    std::array<double, 3> newPosition = (particle.getX()) + deltaT_ * (particle.getV()) + force; // calculate position
                    unsigned int afterCellIndex = translate3DPosTo1D(newPosition);
                    if (translate3DIndTo1D(x,y,z) == afterCellIndex)
                    {
                        particle.setX(newPosition);
                    }
                    else
                    {
                        // the particle has to move to a new cell
                        particle.setX(newPosition);
                        currentCell.erase(currentCell.begin() + i);
                        --i;
                        //indexesToDelete.emplace_back(i);
                    /*    Particle temp = *it;
                        cells.at(beforeCellIndex).erase(it);
                        (*it).setX(newPosition);
                        addBack.push_back(temp);*/
                    }
                }
               /* for (unsigned int i = 0; i < indexesToDelete.size(); ++i) {

                }*/
            }
        }
    }
    // add all particles from addback to the cells
    for (unsigned int i = 0; i < addBack.size(); ++i) {
        cells.at(translate3DPosTo1D(addBack.at(i).getX())).emplace_back(addBack.at(i));
    }
    iterHalo();
}


/*void ParticleContainerLinCel::calculatePosition()
{
    std::cout << "Richtiges calculatePosition\n";
    std::vector<Particle> addBack;
    *//*auto updateLambda = [&](ParticleContainerLinCel::cell::iterator it)
    {
        unsigned int beforeCellIndex = translate3DPosTo1D((*it).getX());
        std::array<double, 3> force = (*it).getF();
        double factor = std::pow(deltaT_, 2) / (2 * (*it).getM());
        force = factor * force;
        std::array<double, 3> newPosition = ((*it).getX()) + deltaT_ * ((*it).getV()) + force; // calculate position
        unsigned int afterCellIndex = translate3DPosTo1D(newPosition);
        if (beforeCellIndex == afterCellIndex)
        {
            (*it).setX(newPosition);
            it++;
        }
        else
        {
            Particle temp = *it;
            cells.at(beforeCellIndex).erase(it);
            (*it).setX(newPosition);
            addBack.push_back(temp);
        }
    };*//*


    iterOverAllParticles(updateLambda);
    for (const auto& pI : addBack)
    {
        unsigned int newIndex = translate3DPosTo1D(pI.getX());
        cells.at(newIndex).push_back(pI);
    }
    iterHalo();
}*/

unsigned int ParticleContainerLinCel::getAmountOfCells() const
{
    return amountOfCells;
}

std::vector<std::vector<Particle>> ParticleContainerLinCel::getCells()
{
    return cells;
}

size_t ParticleContainerLinCel::getAmountOfParticles()
{
    return particles_.size();
}

void ParticleContainerLinCel::calculateForces()
{
    for (auto &p : particles_)
    {
        auto oldForce = p->getF();
        std::array<double, 3> zero = {0.0, 0.0, 0.0};
        p->setF(zero);
        p->setOldF(oldForce);
    }
    iterOverInnerPairs(force_);
    iterBoundary();
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
    if (position[0] >= 0) {
        if (position[0] >= (domainSize_[0] + 1)) {
            xIndex = cellsX - 1;
        }
        xIndex = static_cast<unsigned int>(floor(position[0] / cutoffRadius_)) + 1;
    }
    else {
        xIndex = 0;
    }
    if (position[1] >= 0) {
        if (position[1] >= (domainSize_[1] + 1)) {
            yIndex = cellsY - 1;
        }
        yIndex = static_cast<unsigned int>(floor(position[1] / cutoffRadius_)) + 1;
    }
    else {
        yIndex = 0;
    }
    if (position[2] >= 0) {
        if (position[2] >= (domainSize_[2] + 1)) {
            zIndex = cellsZ - 1;
        }
        zIndex = static_cast<unsigned int>(floor(position[2] / cutoffRadius_)) + 1;
    }
    else {
        zIndex = 0;
    }
    unsigned int index = xIndex + cellsX * yIndex + cellsX * cellsY * zIndex;
    return index;
}

std::array<int, 3> ParticleContainerLinCel::translate1DIndTo3DInd(int index) const
{
    // TODO: overthink this again
    int plane = cellsX * cellsY;
    // calculate the z coordinate
    int z = index / plane;
    index = index % z;
    // calculate the y coordinate
    int y = index / cellsX;
    index % y;
    int x = index;
    return std::array<int, 3>{x, y, z};
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
    auto numberofDimensions = 3;
    return calculateKinEnergy() / getAmountOfParticles() * 2.0 / numberofDimensions;
}

bool ParticleContainerLinCel::affectsForce(int index)
{
    if (index < 0 || index >= 6)
    {
        return false;
    }
    return conditions_.at(index).get()->affectsForce();
}

void ParticleContainerLinCel::iterBoundary2() {

}

std::vector<Particle> ParticleContainerLinCel::getAllParticles() {
    std::vector<Particle> returnVector;
    for (auto &cell : cells) {
        for (int i = 0; i < cell.size(); ++i) {
            returnVector.emplace_back(cell.at(i));
        }

    }
}

void ParticleContainerLinCel::simulateParticles2() {
    auto begin = std::chrono::high_resolution_clock::now();
    this->cells;
    std::vector<Particle> allParticles;// = getAllParticles(); //TODO <-- hieran liegt der Fehler
    for (auto &cell : cells) {
        for (int i = 0; i < cell.size(); ++i) {
            allParticles.emplace_back(cell.at(i));
        }

    }
    while (startTime_ < endTime_)
    {
        /*if (outManager_->outputFiles && iteration_ % outputEveryNIterations_ == 0)
        {*/
        outManager_->plotParticles2(allParticles, iteration_);
        //}
        if (iteration_ % 100 == 0)
        {
            LogManager::infoLog("Iteration {} finished. ({}%)", iteration_, std::round(iteration_ * 10000 / (endTime_ / deltaT_)) / 100);
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
       /* std::cout << "bis hier ok 4\n";

        std::cout << "richtiges simulateParticles\n";*/
        // calculate new x
        calculatePosition();
        // std::cout << "bis hier ok: nach calculatePosition\n";
        // calculate new f
        calculateForces();
        // calculate new v
        calculateVelocity();

        iteration_++;
        startTime_ += deltaT_;
    }
    //  TODO (ADD): Log
    auto end = std::chrono::high_resolution_clock::now();
    size_t diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    std::cout << "Output written, took " + std::to_string(diff) + " milliseconds. (about " + (iteration_ > diff ? std::to_string(iteration_ / diff) + " iter/ms" : std::to_string(diff / iteration_) + " ms/iter") + ") Terminating...\n";

}