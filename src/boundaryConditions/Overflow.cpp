//
// Created by alex on 28.11.23.
//

#include "Overflow.h"

/*
 * direction:
 *  0 == x
 *  1 == y
 *  2 == z
 */
std::function<void(Particle &a)>
Overflow::applyBoundary(std::function<void(Particle &, Particle &)> &forceLambda, double position, int direction) {
    Particle ghostParticle;
    // derive the position of the ghostParticle based on position and direction of the boundary

    // return
}
