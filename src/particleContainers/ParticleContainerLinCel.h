#pragma once
#include "particleContainers/ParticleContainer.h"
#include "boundaryConditions/BoundaryCondition.h"
#include "forces/Force.h"

class ParticleContainerLinCel : public ParticleContainer
{
private:
    using cell = std::vector<Particle>;
    std::vector<cell> cells;
    std::vector<BoundaryCondition> conditions_;
    bool upperModulo = false;
    bool rightModulo = false;
    double amountOfCells = 0;
    /**
     * @brief the cells can be divided into inner, boundary and halo cells
     *
     * the first element of the cells is the cell at the front lower left corner of the domain
     * cells go from:
     * 1. left to right
     * 2. from down to up
     * 3. from front to back
     */
    std::array<double, 3> domainSize_;
    /**
     * 2D cell: cutoffRadius * cutoffRadius * 1
     * 3D cell: cutoffRadius * cutoffRadius * cutoffRadius
     */
    double cutoffRadius_;
    // amount of Cells in each dimension can only be an unsigned integer
    unsigned int cellsX = 0;
    unsigned int cellsY = 0;
    // only needed for the 3D case
    unsigned int cellsZ = 0;

public:

    /**
     * @brief constructor
     *  If the domain size isn't a multiple of the cutoff radius
     *  the last row and/or column of cells will only be a fraction of a cell
     *  for example domainSize{10, 9} with cutoff radius 3
     *  the last column of the cells will be only a (1 * cutoffRadius) cell instead of a (cutoffradius * cutoffRadius) Cell
     * @param domainSize  size of the domain
     * @param bounds  string describing the boundary conditions where each char represents a different side of the domain
     * @param force the force calulation which will be used in this container
     * @param cutoffRadius the cutoff radius according to the Linked Cell algorithm
     */
    ParticleContainerLinCel(double deltaT, double endTime, int writeFrequency, const std::array<double, 3> &domainSize, const std::string &bounds, Force &force, double cutoffRadius);
    /**
     * @brief destructor
     */
    ~ParticleContainerLinCel();

    void iterOverInnerPairs(const std::function<void(Particle &a, Particle &b)> &f) override;
    /**
     * @brief overriding the add method so that particles get added to the correct cell
     * according to their coordinates
     * if a new cell is iterated over, only calculate the forces in the cell itself
     * and the forces between the particles in the current cell and the particles on the right hand side of the cell and
     * over the current cell, also the upperLeft cell has to be checked (a check for the cutoffRadius has to be included)
     * Calc       Calc   Calc
     *      \    |     /
     *       \   |   /
     *        \ Cell -- > Calc
     * @param x_arg position of the particle
	 * @param v_arg velocity of the particle
	 * @param mass mass of the particle
	 * @param type typenumber of the particle
     */
    void add(const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg, double mass, int type) override;

    /**
     * @brief iterate over the particles which are currectly located in the boundary zone and
     * apply the boundary condition, the cells at the corner will be iterated over twice
     * @param None
     * @return void
     */
    void iterBoundary();

    /**
     * @brief iterate over the particles which are currectly located in the halo zone and
     * delete them
     * @param None
     * @return void
     */
    void iterHalo();

    /**
     * @brief overriding the position calculation, so that the particles and cells get updated
     * @param None
     * @return void
     */
    void calculatePosition() override;

    /**
     * @brief overriding the force calculation for the demands of the LinCel container, including iterating over the
     * boundary cells and applying the boundary conditions
     * @param None
     * @return void
     */
    void calculateForces() override;

    /**
     * @brief return the sum of all particles in all the cells of the LinCel container
     * @param None
     * @return void
     */
    unsigned int getAmountOfParticles();

    /**
     * @brief return of how many cells the current LinCel container consists
     * @param None
     * @return void
     */
    unsigned int getAmountOfCells() const;

    /**
     * @brief return a vector of the cells
     * @param None
     * @return void
     */
    std::vector<std::vector<Particle>> getCells();
};
