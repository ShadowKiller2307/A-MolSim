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

    if (domainSize_[2] < 0.0) {
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
    if (fmod(domainSize_[2], cutoffRadius_) != 0.0) {
        depthModulo = true;
    }

    cellsX = static_cast<unsigned int>(ceil(domainSize_[0] / cutoffRadius_) + 2);
    cellsY = static_cast<unsigned int>(ceil(domainSize_[1] / cutoffRadius_) + 2);
    cellsZ = static_cast<unsigned int>(ceil(domainSize_[2]/cutoffRadius_) + 2);
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
    
    for (int x = 1; x < cellsX-1; ++x) {
        for (int y = 1; y < cellsY-1; ++y) {
            for (int z = 1; z < cellsZ-1; ++z) {
                cell &current = cells.at(translate3DIndTo1D(x, y, z));
                for (auto &pi : current) {
                    // right cell
                    if (x != cellsX - 2) {
                        for (auto &pj: cells.at(translate3DIndTo1D(x+1, y, z))) {
                            // check whether pj is within the cutoff radius of pi
                            if (ArrayUtils::L2Norm(pj.getX() - pi.getX()) <= cutoffRadius_) {
                                f(pi, pj);
                            }
                        }
                    }
                    // upper cell
                    if (y != cellsY - 2) {
                        for (auto &pj: cells.at(translate3DIndTo1D(x, y+1, z))) {
                            // check whether pj is within the cutoff radius of pi
                            if (ArrayUtils::L2Norm(pj.getX() - pi.getX()) <= cutoffRadius_) {
                                f(pi, pj);
                            }
                        }
                    }
                    //right upper cell
                    if ((x != cellsX - 2) && (y != cellsY - 2)) {
                        for (auto &pj: cells.at(translate3DIndTo1D(x+1, y+1, z))) {
                            // check whether pj is within the cutoff radius of pi
                            if (ArrayUtils::L2Norm(pj.getX() - pi.getX()) <= cutoffRadius_) {
                                f(pi, pj);
                            }
                        }
                    }
                    // left upper cell
                    if ((x > 1) && (y != cellsY - 2)) {
                        for (auto &pj: cells.at(translate3DIndTo1D(x-1, y+1, z))) {
                            // check whether pj is within the cutoff radius of pi
                            if (ArrayUtils::L2Norm(pj.getX() - pi.getX()) <= cutoffRadius_) {
                                f(pi, pj);
                            }
                        }
                    }
                    //TODO: expand implementation for 3D cases

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

}

void ParticleContainerLinCel::iterHalo()
{
    for (int x = 0; x < cellsX; ++x) {
        for (int y = 0; y < cellsY; ++y) {
            for (int z = 0; z < cellsZ; ++z) {
                if ((x == 0) || (x == cellsX - 1) || (y == 0) || (y == cellsY -1) || (z == 0)
                || (z == cellsZ -1)) {
                    cell &currentCell = cells.at(translate3DIndTo1D(x, y, z));
                    currentCell.clear();
                }
            }
        }
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
