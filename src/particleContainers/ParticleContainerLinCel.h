#pragma once
#include "particleContainers/ParticleContainer.h"
#include "boundaryConditions/BoundaryCondition.h"
#include "forces/Force.h"
#include "boundaryConditions/Outflow.h"
#include "boundaryConditions/Reflecting.h"

enum class BoundaryCondition2
{
    Reflecting2,
    Periodic2,
    Outflow2
};

class ParticleContainerLinCel : public ParticleContainer
{

private:
    using cell = std::vector<Particle>;
    std::vector<cell> cells;
    std::vector<std::unique_ptr<BoundaryCondition>> conditions_;
    std::vector<BoundaryCondition2> conditions2;
    std::vector<Reflecting> reflectingBounds;
    std::vector<Outflow> outflowBounds;
    bool upperModulo = false;
    bool rightModulo = false;
    bool depthModulo = false;
    uint32_t amountOfCells = 0;
    std::vector<std::vector<std::vector<int>>> lookup;
    /**
     * @brief the cells can be divided into inner, boundary and halo cells
     *
     * the first element of the cells is the cell at the front lower left corner of the domain
     * cells go from:
     * 1. left to right
     * 2. from down to up
     * 3. from back to front
     */
    std::array<double, 3> domainSize_;
    /**
     * 2D cell: cutoffRadius * cutoffRadius * 1
     * 3D cell: cutoffRadius * cutoffRadius * cutoffRadius
     */
    double cutoffRadius_;
    // amount of Cells in each dimension can only be an unsigned integer
    uint32_t cellsX = 0;
    uint32_t cellsY = 0;
    // only needed for the 3D case
    uint32_t cellsZ = 0;
    void buildLookUp();
    std::function<void(Particle &)> createPeriodicLambda(int direction, int position);

public:
    bool cellPointerNeedUpdate;

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
    ~ParticleContainerLinCel() = default;

    // void recalculateParticlesinCells();

    void simulateParticles2();

    /**
     * if a new cell is iterated over, only calculate the forces in the cell itself
     * and the forces between the particles in the current cell and the particles on the right hand side of the cell and
     * over the current cell, also the upperLeft cell has to be checked (a check for the cutoffRadius has to be included)
     * Calc       Calc   Calc
     *      \    |     /
     *       \   |   /
     *        \ Cell -- > Calc
     *
     */
    void iterOverInnerPairs(const std::function<void(Particle &a, Particle &b)> &f) override;

    // void iterOverAllParticles(const std::function<void(ParticleContainerLinCel::cell::iterator)> &f);

    void iterOverAllParticles(const std::function<void(ParticleContainerLinCel::cell::iterator)> &f);

    /**
     * @brief overriding the add method so that particles get added to the correct cell
     * according to their coordinates
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

    void iterBoundary2();

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
     * @brief return the number of all particles in all the cells of the LinCel container
     * @param None
     * @return number of particles in the domain
     */
    size_t getAmountOfParticles();

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

    /**
     * @brief translate a 3D cell index to a 1D cell index(for our cells vector)
     * first index (0,0,0) corresponds to the left front corner cell within the domain (so not halo cells)
     * @param x xIndex of the cell
     * @param y yIndex of the cell
     * @param z zIndex of the cell
     * @return 1D index for our cells vector
     */
    unsigned int translate3DIndTo1D(unsigned int x, unsigned int y, unsigned int z) const;

    /**
     * @brief translate 3d coordinates to an index of our cells vector
     * @param position the position as 3D coordinate
     * @return the index in our vector
     */
    unsigned int translate3DPosTo1D(std::array<double, 3> position) const;

    /**
     *
     */
    std::array<int, 3> translate1DIndTo3DInd(int index) const;

    double calculateKinEnergy();

    double calculateTemperature();

    bool affectsForce(int index);

    std::vector<Particle> getAllParticles();
    // std::vector<BoundaryCondition> getBounds();
};
