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
    if (fmod(domainSize_[0], cutoffRadius) != 0.0)
    { // last column of cells are smaller in width
        rightModulo = true;
    }
    if (fmod(domainSize[1], cutoffRadius) != 0.0)
    { // last row of cells are smaller in height
        upperModulo = true;
    }

    cellsX = static_cast<unsigned int>(ceil(domainSize_[0] / cutoffRadius_) + 2);
    cellsY = static_cast<unsigned int>(ceil(domainSize_[1] / cutoffRadius_) + 2);
    // amountOfCells = cellsX * cellsY + 2 * (cellsX + 2) + 2 * cellsY;
    amountOfCells = cellsX * cellsY;
    /**
     *  If the domain size isn't a multiple of the cutoff radius
        the last row and/or column of cells will only be a fraction of a cell
        for example domainSize{10, 9} with cutoff radius 3
        the last column of the cells will be only a (1 * cutoffRadius) cell instead of a (cutoffradius * cutoffRadius) Cell
     */
    for (int i = 0; i < amountOfCells; ++i)
    { // intialize the cells
        cells.emplace_back();
    }
}

void ParticleContainerLinCel::iterOverInnerPairs(const std::function<void(Particle &, Particle &)> &f)
{
    for (int x = 0; x < amountOfCells; ++x)
    {
        if (x < cellsX || x >= (amountOfCells - cellsX) || x % cellsX == 0 || x % cellsX == (cellsX - 1))
        {
            continue; // skip the halo cells
        }
        cell &current = cells.at(x);
        bool rightCell = (x % cellsX) < (cellsX - 2);
        bool upperCell = x < (amountOfCells - 2 * cellsX);
        bool rightUpperCell = rightCell and upperCell;
        bool leftUpperCell = ((x % cellsX) >= 2) and upperCell;
        /*
         * store the right and upper cells in a reference, if no right or upper cell exists
         * store the current cell in the variable(so that the variables always have a value)
         */
        cell &currentRight = rightCell ? cells.at(x + 1) : cells.at(x);
        cell &currentUpper = upperCell ? cells.at(x + cellsX) : cells.at(x);
        cell &currentUpperRight = rightUpperCell ? cells.at(x + cellsX + 1) : cells.at(x);
        cell &currentUpperLeft = leftUpperCell ? cells.at(x + cellsX - 1) : cells.at(x);
        for (int i = 0; i < current.size(); ++i)
        {
            Particle &pi = current.at(i);
            // forces within the current cell
            for (int j = i + 1; j < current.size(); ++j)
            {
                Particle &pj = current.at(j);
                // every particle in the same cell is within the cutoff radius
                f(pi, pj);
            }
            if (rightCell)
            {
                for (auto &pj : currentRight)
                {
                    // check whether pj is within the cutoff radius of pi
                    if (ArrayUtils::L2Norm(pj.getX() - pi.getX()) <= cutoffRadius_)
                    {
                        f(pi, pj);
                    }
                }
            }
            if (upperCell)
            {
                for (auto &pj : currentUpper)
                {
                    // check whether pj is within the cutoff radius of pi
                    if (ArrayUtils::L2Norm(pj.getX() - pi.getX()) <= cutoffRadius_)
                    {
                        f(pi, pj);
                    }
                }
            }
            if (rightUpperCell)
            {
                for (auto &pj : currentUpperRight)
                {
                    // check whether pj is within the cutoff radius of pi
                    if (ArrayUtils::L2Norm(pj.getX() - pi.getX()) <= cutoffRadius_)
                    {
                        f(pi, pj);
                    }
                }
            }
            if (leftUpperCell)
            {
                for (auto &pj : currentUpperLeft)
                {
                    // check whether pj is within the cutoff radius of pi
                    if (ArrayUtils::L2Norm(pj.getX() - pi.getX()) <= cutoffRadius_)
                    {
                        f(pi, pj);
                    }
                }
            }
        }
    }
}

void ParticleContainerLinCel::iterOverAllParticles(const std::function<void(Particle &)> &f)
{
    for (size_t x = 1; x < cellsX - 1; x++)
    {
        for (size_t y = 1; y < cellsY - 1; y++)
        {
            for (size_t z = 1; z < cellsZ - 1; z++)
            {
                for (auto &p : cells[translate3DIndTo1D(x, y, z)])
                {
                    f(p);
                }
            }
        }
    }
}

void ParticleContainerLinCel::iterBoundary()
{
    unsigned int i = cellsX + 1;
    if (conditions_[2].affectsForce())
    {
        // from left to right
        for (i; i < cellsX - 2; ++i)
        { // lower row
            cell &currentCell = cells.at(i);
            for (auto &pi : currentCell)
            {
                conditions_[2].applyBoundCondition(pi);
            }
        }
    }
    if (conditions_[1].affectsForce())
    {
        // from down to up
        for (int j = 0; j < cellsY - 2; ++j)
        {
            cell &currentCell = cells.at(i);
            for (auto &pi : currentCell)
            {
                conditions_[1].applyBoundCondition(pi);
            }
            if (rightModulo)
            {
                cell &additionalCell = cells.at(i - 1);
                for (auto &pi : additionalCell)
                {
                    if ((domainSize_[0] - pi.getX()[0]) < cutoffRadius_)
                    {
                        conditions_[1].applyBoundCondition(pi);
                    }
                }
            }
            i += cellsX;
        }
    }
    i -= cellsX;

    if (conditions_[3].affectsForce())
    {
        // from right to left
        for (int j = 0; j < cellsX - 2; ++j)
        {
            cell &currentCell = cells.at(i);
            for (auto &pi : currentCell)
            {
                conditions_[3].applyBoundCondition(pi);
            }
            if (upperModulo)
            {
                cell &additionalCell = cells.at(i - cellsX);
                for (auto &pi : additionalCell)
                {
                    if ((domainSize_[1] - pi.getX()[1]) < cutoffRadius_)
                    {
                        conditions_[3].applyBoundCondition(pi);
                    }
                }
            }
            --i;
        }
    }
    ++i;
    if (conditions_[0].affectsForce())
    {
        // from up to down
        for (int j = 0; j < cellsY - 2; ++j)
        {
            cell &currentCell = cells.at(i);
            for (auto &pi : currentCell)
            {
                conditions_[0].applyBoundCondition(pi);
            }
            i -= cellsX;
        }
    }
}

void ParticleContainerLinCel::iterHalo()
{
    int i = 0;
    // from left to right
    for (i; i < cellsX; ++i)
    {
        cell &currentCell = cells.at(i);
        currentCell.clear(); // delete all particles from the current halo cell
    }
    --i;
    // from down to up
    for (int j = 0; j < cellsY - 1; ++j)
    {
        i += cellsX;
        cell &currentCell = cells.at(i);
        currentCell.clear();
    }
    // from right to left
    for (int j = 0; j < cellsX - 1; ++j)
    {
        --i;
        cell &currentCell = cells.at(i);
        currentCell.clear();
    }
    // from up to down
    for (int j = 0; j < cellsY - 2; ++j)
    {
        i = i - cellsX;
        cell &currentCell = cells.at(i);
        currentCell.clear();
    }
}

void ParticleContainerLinCel::add(const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg, double mass,
                                  int type)
{
    // compute the cell to which the particle will be added
    if (x_arg[0] <= domainSize_[0] || x_arg[1] <= domainSize_[1])
    {
        // cells.at(translate3DPosTo1D(x_arg)).emplace_back(x_arg, v_arg, mass, type);
        double xIndex = trunc(x_arg[0] / cutoffRadius_);
        double yIndex = trunc(x_arg[1] / cutoffRadius_);
        double index = static_cast<unsigned int>(xIndex + cellsX * yIndex) + 1 + cellsX; // The "+1 + cellsX" is added because of the halo cells
        cells.at(index).emplace_back(x_arg, v_arg, mass, type);
    }
}

// TODO: check calculatePosition for errors
void ParticleContainerLinCel::calculatePosition()
{
    int i = 0;
    for (int j = 0; j < cells.size(); ++j) // loop um durch alle cells durchzugehen
    {
        auto &currentCell = cells.at(j); // aktuelle cell
        for (int x = 0; x < currentCell.size(); ++x)
        {
            // for (auto &p : currentCell) {
            auto &p = currentCell.at(x); // aktuelles particle in der aktuellen cell
            std::array<double, 3> force = p.getF();
            double factor = std::pow(deltaT_, 2) / (2 * p.getM());
            force = factor * force;
            std::array<double, 3> newPosition = p.getX() + deltaT_ * p.getV() + force; // berechne die position
            // check whether the particle left the current cell
            double xIndexNewCell; // soll den x-index der cell beschreiben in der das partikel mit der neuen position sein wird
            double yIndexNewCell;
            if (newPosition[0] < 0.0)
            {
                xIndexNewCell = 0;
            }
            else
            {
                xIndexNewCell = trunc(newPosition[0] / cutoffRadius_) + 1; // loscht die Nachkommastellen
            }
            if (newPosition[1] < 0.0)
            {
                yIndexNewCell = 0;
            }
            else
            {
                yIndexNewCell = trunc(newPosition[1] / cutoffRadius_) + 1; // loscht die Nachkommastellen
            }
            double indexNewCell = xIndexNewCell + cellsX * yIndexNewCell;
            if (indexNewCell >= amountOfCells || indexNewCell < 0)
            {
                currentCell.erase(currentCell.begin() + x);
            }
            else if (j != indexNewCell)
            { // Particle has to be deleted in the old cell and added to the new cell
                currentCell.erase(currentCell.begin() + x);
                p.setX(newPosition);
                cells.at(indexNewCell).emplace_back(p);
            }
            i++;
        }
    }
}

unsigned int ParticleContainerLinCel::getAmountOfCells() const
{
    return amountOfCells;
}

std::vector<std::vector<Particle>> ParticleContainerLinCel::getCells()
{
    return cells;
}

unsigned int ParticleContainerLinCel::getAmountOfParticles()
{
    unsigned int returnValue = 0;
    for (auto cell : cells)
    {
        returnValue += cell.size();
    }
    return returnValue;
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
    iterBoundary();
}

int ParticleContainerLinCel::translate3DIndTo1D(int x, int y, int z)
{
    int index = (x + 1) + cellsX * (y + 1) + cellsX * cellsY * (z + 1);
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

unsigned int ParticleContainerLinCel::translate3DPosTo1D(std::array<double, 3> position)
{
    unsigned int xIndex = static_cast<unsigned int>(floor(position[0] / cutoffRadius_)) + 1;
    unsigned int yIndex = static_cast<unsigned int>(floor(position[1] / cutoffRadius_)) + 1;
    unsigned int zIndex = static_cast<unsigned int>(floor(position[2] / cutoffRadius_)) + 1;
    unsigned int index = xIndex + cellsX * yIndex + cellsX * cellsY * zIndex;
    return index;
}

std::array<int, 3> ParticleContainerLinCel::translate1DIndTo3DInd(int index)
{
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
    auto energyLambda = [&energy](Particle &p)
    {
        energy += p.getM() * std::inner_product(p.getV().begin(), p.getV().end(), p.getV(), 0.0) / 2;
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
