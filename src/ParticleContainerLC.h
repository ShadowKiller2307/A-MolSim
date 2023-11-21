/**
 * @brief This class implements the ParticleContainer using the linked cell algortihm
 */

#pragma once
#include "ParticleContainer.h"
#include "Particle.h"
#include "cmath"

using cell = std::vector<Particle>;

class ParticleContainerLC : ParticleContainer {
private:
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
    double cellsX = domainSize[0]/cutoffRadius;
    double cellsY = domainSize[1]/cutoffRadius;
    double cellsZ = domainSize[2]/cutoffRadius;

public:
    ParticleContainerLC(std::array<double, 3> domainSize, double cutoffRadius);

    void iterOverPairs(const std::function<void(Particle &a, Particle &b)> &f);
    /**
     * @brief iterate over the particles which are currectly located in the boundary zone
     */
    void iterBoundary();
    /**
     * @brief iterate over the particles which are currectly located in the halo zone
     */
     void iterHalo();

//
//     void find_index(double coordinate) {
//         double counter = 0.0;
//         for (int index = 0; index < ; ++i) {
//             if (counter + cutoffRadius > coordinate)
//         }
//
//     }
     /**
      * @brief method to add an Particle to an LC container, thus adding it to the according cell
      */
      void add(Particle &a) {
          // compute the cell to which the particle will be added
          double xIndex = trunc(a.getX()[0] / cutoffRadius);
          double yIndex = trunc(a.getX()[1]/cutoffRadius);
          double zIndex = trunc(a.getX()[2]/cutoffRadius);
          double index = xIndex + yIndex * cellsY; //TODO complete the index calculation

      }
};
