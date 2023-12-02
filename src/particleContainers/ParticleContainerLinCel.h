#pragma once
#include "particleContainers/ParticleContainer.h"
#include "boundaryConditions/BoundaryCondition.h"
using cell = std::vector<Particle>;

class ParticleContainerLinCel : public ParticleContainer
{
private:
    using cell = std::vector<Particle>;
	std::vector<cell> cells;
    std::vector<BoundaryCondition> conditions_;
    bool upperModulo = false;
    bool rightModulo = false;

    double amountOfCells = 0;
    // verschiebe die Rechnung in den Konstruktor
    /*
    = floor(domainSize_[0]/cutoffRadius_) * floor(domainSize_[1]/cutoffRadius_) // (domainSize_[2]/cutoffRadius_) //inner+boundary cells
                           + 2 * cellsX + 2 * (cellsY-2); // adding the halo cells*/
    /**
    * @brief the cells can be divided into inner, boundary and halo cells
    *
    * the first element of the cells is the cell at the front lower left corner of the domain
    * cells go from:
    * ->left to right
    * ->front to back
    * ->from down to up
    */
    std::array<double, 3> domainSize_;
    /**
     * 2D cell: cutoffRadius * cutoffRadius * 1
     * 3D cell: cutoffRadius * cutoffRadius * cutoffRadius
     */
    double cutoffRadius_;
    /*
     * amount of Cells in each dimension can only be an unsigned integer
     */
    unsigned int cellsX = 0;
    unsigned int cellsY = 0;
    // only needed for the 3D case
    unsigned int cellsZ = 0;

public:
	ParticleContainerLinCel(double deltaT, double endTime, std::array<double, 3> domainSize, double cutoffRadius, std::vector<BoundaryCondition> &conditions);


	~ParticleContainerLinCel();


    void iterOverInnerPairs(const std::function<void(Particle &a, Particle &b)> &f) override;
    /**
     * @brief iterate over the particles which are currectly located in the boundary zone
     *  @param boundaryLambda an array of 4 BoundaryConditions so that a different boundary condition can be applied to each side
     */
    void add(const std::array<double, 3> &x_arg, const std::array<double, 3> &v_arg, double mass, int type) override;

    /**
     * @brief iterate over the particles which are currectly located in the halo zone
     */
    void iterBoundary(std::array <const std::function<void(Particle &, Particle &)>, 4> &boundaryLambda);

    /**
     * @brief iterate over the particles which are currectly located in the halo zone
     */
    void iterHalo();

    void calculatePosition() override;

    unsigned int getAmountOfCells();

    std::vector<std::vector<Particle>> getCells();
};
