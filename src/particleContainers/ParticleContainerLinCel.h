#pragma once
#include "particleContainers/ParticleContainer.h"
#include "forces/Force.h"
#include "Thermostat.h"
#include <array>

enum class BoundaryCondition
{
    Reflecting,
    Periodic,
    Outflow
};

class ParticleContainerLinCel : public ParticleContainer
{

private:
    /// @brief vector holding all particles belonging to that cell
    using cell = std::vector<Particle>;
    /// @brief 1D vector of cells
    std::vector<cell> cells;
    /**
     * @brief the boundary conditions of each face of the domain
     * 0: towards negative x
     * 1: towards positive x
     * 2: towards negative y
     * 3: towards positive y
     * 4: towards negative z
     * 5: towards positive z
     */
    std::vector<BoundaryCondition> conditions;
    /// @brief how many cells there are, equivalent to cells.size()
    uint32_t amountOfCells = 0;
    /// @brief currently unused lookuptable to speed up 3D to 1D translations
    std::vector<std::vector<std::vector<int>>> lookup;
    unsigned int mUpdates;
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
    /// @brief how many cells there are in each dimension respactively
    uint32_t cellsX, cellsY, cellsZ = 0;
    /// @brief a way to access cellsX, cellsY, cellsZ by index, 0 is for X, 1 for Y, 2 for Z
    std::array<uint32_t, 3> cellCountByIndex;

    /// @brief specifies if the simulation uses a thermostat
    bool useThermostat;
    /// @brief after how many iterations the thermostat should be applied each time
    unsigned int nThermostat = 100;
    /// @brief specifies if the velocity scaling is gradual or direct
    bool isGradual;
    // the thermostat for this container
    // Thermostat thermostat;
    double tempTarget;
    double maxDiff;
    void buildLookUp();
    double initT = 0;
    double gGrav;

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
    ParticleContainerLinCel(double deltaT, double endTime, int writeFrequency,
                            const std::array<double, 3> &domainSize,
                            const std::string &bounds, double cutoffRadius,
                            bool useThermostat = false, unsigned int nThermostat = 100,
                            bool isGradual = true, double initT = 0,
                            double tempTarget = 20,
                            double maxDiff = 0.5,
                            double gGrav = 0);

    /**
     * @brief destructor for Linked Cells Container
     */
    ~ParticleContainerLinCel() = default;

    /**
     * @brief overriding the add method so that particles get added to the correct cell
     * according to their coordinates
     * @param x_arg position of the particle
     * @param v_arg velocity of the particle
     * @param mass mass of the particle
     * @param type typenumber of the particle
     */
    void add(const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg, double mass, int type, double epsilon = 5.0, double sigma = 1.0) override;

    /**
     * @brief adds a fully constructed particle to the container
     *
     * @param p the particle to add
     * @return void
     */
    void addCompleteParticle(Particle &p) override;

    /// @brief runs the simulation loop
    void simulateParticles() override;

    /**
     * @brief overriding the force calculation for the demands of the LinCel container, including iterating over the
     * boundary cells and applying the boundary conditions
     * @param None
     * @return void
     */
    void calculateForces() override;

    /**
     * @brief calculates the forces between the particles of two cells which are passed as coordinates.
     * If the cells are not next to each other the particls will be temporarily be moved in order to get the correct forces and use the cutoffRadius
     *
     * @param myCoordinates coordinates of the first cell
     * @param otherCoordinates coordinates of the second cell
     * @return void
     */
    void calculateForcesWithIndices(std::array<uint32_t, 3> &myCoordinates, std::array<uint32_t, 3> &otherCoordinates);

    /**
     * @brief updates the positions of all particles according to their forces and velocities
     * @param None
     * @return void
     */
    void calculatePosition() override;

    /**
     * @brief updates the velocities of all particles according to their forces
     * @param None
     * @return void
     */
    void calculateVelocity() override;

    /**
     * @brief iterates over all particles individually and applies the function to each one. The function itself is responsible for advancing the iterator or not
     * @param f the function to be applied
     * @return void
     */
    void iterOverAllParticles(const std::function<void(ParticleContainerLinCel::cell::iterator)> &f);

    /**
     * @brief if a new cell is iterated over, only calculate the forces in the cell itself
     * and the forces between the particles in the current cell and the particles on the right hand side of the cell and
     * over the current cell, also the upperLeft cell has to be checked (a check for the cutoffRadius has to be included)
     * Calc       Calc   Calc
     *      \    |     /
     *       \   |   /
     *        \ Cell -- > Calc
     *
     */
    void iterOverInnerPairs(const std::function<void(Particle &a, Particle &b)> &f) override;

    std::function<void(uint32_t x, uint32_t y, uint32_t z)> createReflectingLambdaBoundary(int direction, int position);
    std::function<void(uint32_t x, uint32_t y, uint32_t z)> createPeriodicLambdaBoundary();
    /**
     * @brief iterate over the particles which are located in the boundary zone and
     * apply the boundary condition, the cells at the corner will be iterated over twice
     * @param None
     * @return void
     */
    void iterBoundary();

    void iterBoundary2();

    std::function<void(uint32_t x, uint32_t y, uint32_t z)> createOutflowLambdaHalo();
    std::function<void(uint32_t x, uint32_t y, uint32_t z)> createPeriodicLambdaHalo();
    std::function<void(uint32_t x, uint32_t y, uint32_t z)> createPeriodicLambdaHalo2();

    /**
     * @brief iterate over the particles which are located in the halo zone and
     * delete or move them according to the boundary condition
     * @param None
     * @return void
     */
    void iterHalo();

    void iterHalo2();

    /**
     * @brief translate a 3D cell index to a 1D cell index(for our cells vector)
     * first index (0,0,0) corresponds to the left front corner cell within the domain (so not halo cells)
     * @param x xIndex of the cell
     * @param y yIndex of the cell
     * @param z zIndex of the cell
     * @return 1D index for the cells vector
     */
    unsigned int translate3DIndTo1D(unsigned int x, unsigned int y, unsigned int z) const;

    /**
     * @brief translate 3d coordinates to an index of our cells vector
     * @param position the position as 3D coordinate
     * @return the index in the vector
     */
    unsigned int translate3DPosTo1D(std::array<double, 3> position) const;

    /**
     * @brief return the number of all particles in all the cells of the LinCel container
     * @param None
     * @return number of particles in the domain
     */
    size_t getAmountOfParticles() const override;

    /**
     * @brief return all particles stored in a single vector
     * @param None
     * @return std::vector<Particle>
     */
    std::vector<Particle> getParticles();

    /**
     * @brief return the number of cells the current LinCel container consists of
     * @param None
     * @return number of cells as a uint32_t
     */
    unsigned int getAmountOfCells() const;

    /**
     * @brief return the cells as a vector of vectors
     * @param None
     * @return std::vector<std::vector<Particle>>& cells
     */
    std::vector<std::vector<Particle>> &getCells();

    /**
     * @brief return the boundary conditions as a vector
     * @param None
     * @return const std::vector<BoundaryCondition>& boundary conditions
     */
    const std::vector<BoundaryCondition> &getConditions();

    const std::array<double, 3> &getDomainSize();

    const double getCutOffRadius();

    void addGravitationalForce();

    /////////////////////////////THERMOSTAT BEGIN //////////////////////////////////////////////////////////////////////

    double calculateKinEnergy();

    double calculateTemperature();

    void scaleVelocity(double currentTemp, double newTemp);

    //////////////////////////////THERMOSTAT END //////////////////////////////////////////////////////////////////////
};
