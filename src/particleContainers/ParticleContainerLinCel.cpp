#include "particleContainers/ParticleContainerLinCel.h"
#include "boundaryConditions/Reflecting.h"
#include "boundaryConditions/Outflow.h"
#include "logOutputManager/LogManager.h"
#include <iostream>

/**
 * "rrrrrr"
 * "oooooo"
 * "pppppp"
 */

ParticleContainerLinCel::ParticleContainerLinCel(double deltaT, double endTime, int writeFrequency, const std::array<double, 3> &domainSize, const std::string &bounds, Force &force, double cutoffRadius) : ParticleContainer(deltaT, endTime, writeFrequency, force.innerPairs())
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
        const double pos = i % 2 == 0 ? 0.0 : domainSize[i / 2];
        const int dir = i / 2;
        auto f = force.boundaryPairs();
        if (c == 'r')
        {
            conditions_.push_back(Reflecting(pos, dir, f));
        }
        else if (c == 'o')
        {
            conditions_.push_back(Outflow(pos, dir, f));
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
    for (int i = 0; i < amountOfCells; ++i)
    { // initialize the cells
        cells.emplace_back();
    }
}

void ParticleContainerLinCel::iterOverInnerPairs(const std::function<void(Particle &, Particle &)> &f)
{

    for (int x = 1; x < cellsX - 1; ++x)
    {
        for (int y = 1; y < cellsY - 1; ++y)
        {
            for (int z = 1; z < cellsZ - 1; ++z)
            {
                cell &current = cells.at(translate3DIndTo1D(x, y, z));
                for (auto &pi : current)
                {
                    // right cell
                    if (x != cellsX - 2)
                    {
                        for (auto &pj : cells.at(translate3DIndTo1D(x + 1, y, z)))
                        {
                            // check whether pj is within the cutoff radius of pi
                            if (ArrayUtils::L2Norm(pj->getX() - pi->getX()) <= cutoffRadius_)
                            {
                                f(*pi, *pj);
                            }
                        }
                    }
                    // upper cell
                    if (y != cellsY - 2)
                    {
                        for (auto &pj : cells.at(translate3DIndTo1D(x, y + 1, z)))
                        {
                            // check whether pj is within the cutoff radius of pi
                            if (ArrayUtils::L2Norm(pj->getX() - pi->getX()) <= cutoffRadius_)
                            {
                                f(*pi, *pj);
                            }
                        }
                    }
                    // right upper cell
                    if ((x != cellsX - 2) && (y != cellsY - 2))
                    {
                        for (auto &pj : cells.at(translate3DIndTo1D(x + 1, y + 1, z)))
                        {
                            // check whether pj is within the cutoff radius of pi
                            if (ArrayUtils::L2Norm(pj->getX() - pi->getX()) <= cutoffRadius_)
                            {
                                f(*pi, *pj);
                            }
                        }
                    }
                    // left upper cell
                    if ((x > 1) && (y != cellsY - 2))
                    {
                        for (auto &pj : cells.at(translate3DIndTo1D(x - 1, y + 1, z)))
                        {
                            // check whether pj is within the cutoff radius of pi
                            if (ArrayUtils::L2Norm(pj->getX() - pi->getX()) <= cutoffRadius_)
                            {
                                f(*pi, *pj);
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
}

void ParticleContainerLinCel::iterOverAllParticles(const std::function<void(Particle &, size_t)> &f)
{
    for (size_t x = 1; x < cellsX - 1; x++)
    {
        for (size_t y = 1; y < cellsY - 1; y++)
        {
            for (size_t z = 1; z < cellsZ - 1; z++)
            {
                cell &c = cells[translate3DIndTo1D(x, y, z)];
                for (size_t i = 0; i < c.size(); ++i)
                {
                    f(*c[i], i);
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
 *
 */
void ParticleContainerLinCel::iterBoundary()
{
    for (int z = 1; z < cellsZ - 1; ++z)
    {
        if (z == 1 || z == cellsZ - 2)
        { // all cells in this plane are boundary cells
        }
        else
        { // iterate over the outer ring of inner cells for each inner z index
            // from down left corner to down right corner
            if (conditions_[2].affectsForce())
            {
                for (int x = 1; x < cellsX - 1; ++x)
                {
                    cell &currentCell = cells.at(translate3DIndTo1D(x, 1, z));
                    for (auto &pi : currentCell)
                    {
                        conditions_[2].applyBoundCondition(*pi);
                    }
                }
            }
            // from down right corner to top right corner
            if (conditions_[1].affectsForce())
            {
                for (int y = 1; y < cellsY - 1; ++y)
                {
                    cell &currentCell = cells.at(translate3DIndTo1D(cellsX - 2, y, z));
                    for (auto &pi : currentCell)
                    {
                        conditions_[1].applyBoundCondition(*pi);
                    }
                }
            }
            if (conditions_[3].affectsForce())
            {
                // from top right corner to top left corner
                for (int x = cellsX - 2; x > 0; --x)
                {
                    cell &currentCell = cells.at(translate3DIndTo1D(x, cellsY - 2, z));
                    for (auto &pi : currentCell)
                    {
                        conditions_[3].applyBoundCondition(*pi);
                    }
                }
            }
            // from top left corner to down left corner
            if (conditions_[0].affectsForce())
            {
                for (int y = cellsY - 2; y > 0; --y)
                {
                    cell &currentCell = cells.at(translate3DIndTo1D(1, y, z));
                    for (auto &pi : currentCell)
                    {
                        conditions_[0].applyBoundCondition(*pi);
                    }
                }
            }
        }
    }
}

void ParticleContainerLinCel::iterHalo()
{
    for (int x = 0; x < cellsX; ++x)
    {
        for (int y = 0; y < cellsY; ++y)
        {
            for (int z = 0; z < cellsZ; ++z)
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
        // add the particle
        Particle *p = &particles_.emplace_back(x_arg, v_arg, mass, type);
        // compute the cell to which the particle will be added
        cells.at(translate3DPosTo1D(x_arg)).emplace_back(p);
    }
}

void ParticleContainerLinCel::calculatePosition()
{
    std::vector<Particle *> addBack;
    auto updateLambda = [&](Particle &p, size_t i)
    {
        int beforeCellIndex = translate3DPosTo1D(p.getX());
        std::array<double, 3> force = p.getF();
        double factor = std::pow(deltaT_, 2) / (2 * p.getM());
        force = factor * force;
        std::array<double, 3> newPosition = p.getX() + deltaT_ * p.getV() + force; // calculate position
        int afterCellIndex = translate3DPosTo1D(newPosition);
        if (beforeCellIndex == afterCellIndex)
        {
            p.setX(newPosition);
        }
        else
        {
            cells[beforeCellIndex].erase(cells[beforeCellIndex].begin() + i);
            p.setX(newPosition);
            addBack.push_back(&p);
        }
    };
    iterOverAllParticles(updateLambda);
    for (auto &p : addBack)
    {
        int newIndex = translate3DPosTo1D(p->getX());
        cells.at(newIndex).push_back(p);
    }
}

unsigned int ParticleContainerLinCel::getAmountOfCells() const
{
    return amountOfCells;
}

std::vector<std::vector<Particle *>> ParticleContainerLinCel::getCells()
{
    return cells;
}

unsigned int ParticleContainerLinCel::getAmountOfParticles()
{
    return particles_.size();
}

void ParticleContainerLinCel::calculateForces()
{
    for (auto &cell : cells)
    {
        for (auto &p : cell)
        {
            auto oldForce = p->getF();
            std::array<double, 3> zero = {0.0, 0.0, 0.0};
            p->setF(zero);
            p->setOldF(oldForce);
        }
    }
    iterOverInnerPairs(force_);
    iterBoundary();
}

int ParticleContainerLinCel::translate3DIndTo1D(int x, int y, int z)
{
    int index = x + cellsX * y + cellsX * cellsY * z;
    /*
     * LGS: but sadly not enough equations
     * index-1 = x + cellsX*(y+1) + cellsX*cellsY*(z+1)
     * index-1 = x + cellsX(y+1+cellsY+z+1)
     * index-1 = x + cellsX*y + cellsX + cellsX*cellsY + cellsX*z + cellsX
     * index-1 - cellsX - cellsX*cellsY - cellsX = x+ cellsX*y + cellsX*z
     *
     */
    return index;
}

unsigned int ParticleContainerLinCel::translate3DPosTo1D(std::array<double, 3> position) const
{
    unsigned int xIndex = static_cast<unsigned int>(floor(position[0] / cutoffRadius_)) + 1;
    unsigned int yIndex = static_cast<unsigned int>(floor(position[1] / cutoffRadius_)) + 1;
    unsigned int zIndex = static_cast<unsigned int>(floor(position[2] / cutoffRadius_)) + 1;
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
    auto energyLambda = [&energy](Particle &p, size_t i)
    {
        energy += p.getM() * std::inner_product(p.getV().begin(), p.getV().end(), p.getV().begin(), 0.0) / 2;
    };
    iterOverAllParticles(energyLambda);
    return energy;
}

double ParticleContainerLinCel::calculateTemperature()
{
    auto numberofDimensions = 3;
    return calculateKinEnergy() / getAmountOfParticles() * 2.0 / numberofDimensions;
}

ParticleContainerLinCel::~ParticleContainerLinCel() = default;
