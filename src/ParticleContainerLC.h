/**
 * @brief This class implements the ParticleContainer using the linked cell algortihm
 */

#pragma once
#include "ParticleContainer.h"
#include "Particle.h"
#include "cmath"

class ParticleContainerLC : public ParticleContainer {

public:
    ParticleContainerLC(std::array<double, 3> domainSize, double cutoffRadius);

    /**
     * @brief iterate over the particles which are currectly located in the boundary zone
     *  @param boundaryLambda an array of 4 BoundaryConditions so that a different boundary condition can be applied to each side
     */
    void iterBoundary(std::array <const std::function<void(Particle &, Particle &)>, 4> &boundaryLambda);
    /**
     * @brief iterate over the particles which are currectly located in the halo zone
     */
     void iterHalo();

    void add(Particle &a) override;

    void iterOverPairs(const std::function<void(Particle &a, Particle &b)> &forceLambda) override;
    //void iterOverPairs2D(const std::function<void(Particle &a, Particle &b)> &forceLambda) override;
    void calculatePosition() override;
    void calculateVelocity() override;
    void calculateCellIndex(double xPos, double yPos, double zPos);

  //  void getSize() override;
};
