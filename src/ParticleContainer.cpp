#include "ParticleContainer.h"
#include "Particle.h"
#include "ForceV1.h"
#include "LennardJonesForce.h"


void ParticleContainer::setDeltaTTwo(double deltaT) {
    this->deltaTTwo = deltaT;
}

double ParticleContainer::getDeltaTwo() {
    return this->deltaTTwo;
}

/*
 * The implementation is needed in order to build the project, thus an empty definition whereas the
 * actual definitions are implemented in the subclasses
 */
void ParticleContainer::iterOverPairs(const std::function<void(Particle &, Particle &)> &f) {}

void ParticleContainer::add(Particle &a) {}

// https://sourcemaking.com/design_patterns/strategy/cpp/1 Looked here how the strategy pattern works
void ParticleContainer::setForceCalculator(int mode)
{

    ParticleContainer::debugLog("Setting the mode for the forces to {}.", mode);
    if (mode == 0)
    {
        forceCalculator = new ForceV1();
    }
    else if (mode == 1)
    {
        forceCalculator = new LennardJonesForce{5, 1};
    }
    else
    {
        forceCalculator = new ForceV1(); // If the input mode isn't defined, the forceCalculator is set to ForceV1
    }
}

void ParticleContainer::calculateVelocity() {

}

void ParticleContainer::calculatePosition() {

}
