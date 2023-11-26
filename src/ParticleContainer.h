#pragma once

#include "LogManager.h"
#include "spdlog/spdlog.h"
#include "Particle.h"

class ForceCalculator; //needs to be defined so that the compiler doesn't throw an error

class ParticleContainer
{
protected:
    ForceCalculator *forceCalculator;

    double deltaTTwo;
    std::vector<Particle> particles;
    using cell = std::vector<Particle>;
    double amountOfCells = (domainSize[0]/cutoffRadius) * (domainSize[1]/cutoffRadius) * (domainSize[2]/cutoffRadius)
                             + 2 * cellsX + 2 * (cellsY-2); // adding the halo cells
    /**
    * @brief the cells can be divided into inner, boundary and halo cells
    *
    * the first element of the cells is the cell at the front lower left corner of the domain
    * cells go from:
    * ->left to right
    * ->front to back
    * ->from down to up
    */
    std::vector<cell> cells;
    std::array<double, 3> domainSize;
    /**
     * 2D cell: cutoffRadius * cutoffRadius * 1
     * 3D cell: cutoffRadius * cutoffRadius * cutoffRadius
     */
    double cutoffRadius;
    /*
     * amount of Cells in each dimension can only be integer
     */
    int cellsX = static_cast<int> (domainSize[0]/cutoffRadius);
    int cellsY = static_cast<int> (domainSize[1]/cutoffRadius);
    int cellsZ = static_cast<int> (domainSize[2]/cutoffRadius);
    //TODO: maybe use a std::variant to represent the different particles implementations

public:



    /**
    * @brief sets the deltaT for the container
    * @param deltaT the deltaT passed by the user or the default value
    * @return void
    */
    void setDeltaTTwo(double deltaT);


    /**
     * @brief calculate the velocity for all particles
     * @param None
     * @return void
     */
    virtual void calculateVelocity();

    /**
     * @brief calculate the position for all particles
     * @param None
     * @return void
     */
    virtual void calculatePosition();

    void setForceCalculator(int mode);

    double getDeltaTwo();

    /**
     * @brief iterate over every particle pair in the container and apply the lambda function
     * @param f the lambda function
     */
    //virtual void iterOverPairs3D(const std::function<void(Particle &a, Particle &b)> &f);
    //virtual void iterOverPairs2D(const std::function<void(Particle &a, Particle &b)> &f);
    virtual void iterOverPairs(const std::function<void(Particle &a, Particle &b)> &f);
    /**
     * @brief method to add an Particle to an LC container, thus adding it to the according cell
     */
    virtual void add(Particle &a);

    std::vector<Particle> &getParticles();


    /**
     * @brief sets the particles for the container
     * @param particles1 the particle vector to be set for the container
     * @return void
     */
    void setParticles(const std::vector<Particle> &particles1);

    /**
     * @brief calculate the force for all particles
     * @param None
     * @return void
     */
    void calculateForces();

    /**
    * @brief method to return the overall size of the container
    */
   // virtual void getSize();
};



