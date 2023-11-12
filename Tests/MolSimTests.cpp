/**
 * @brief Class for Unit Tests which will check the correctness of our code
 */

#include <gtest/gtest.h>
#include "../src/ParticleContainer.h"
#include "../src/HelperFunctions.h"
#include "../src/Particle.h"
#include "../src/ParticleGenerator.h"
#include "../src/ForceV1.h"
#include "../src/LennardJonesForce.h"

// Difference ASSERT vs EXPECT macros
// ASSERT -> fatal failures
// EXPECT -> nonfatal failures


// unit Test for the particle container
TEST(ParticleContainerTest, testGetParticles) {
    ParticleContainer container;
    std::vector particles{Particle(0), Particle(1)};
    container.setParticles(particles);
    EXPECT_EQ(2, container.getParticles()->size());
}

//void instantiateCuboid(ParticleContainer& container, std::array<double, 3> llfc, std::array<unsigned int, 3> particlePerDimension, double h, double mass,
//                       std::array<double, 3> particleVelocity);
TEST(ParticleGeneratorTest, testGenerateParticlesGenerator) {
    // Instantiate a generator and container for the instantiateCuboid function
    ParticleGenerator particleGenerator;
    ParticleContainer container;
    particleGenerator.instantiateCuboid(container, {0.0, 0.0, 0.0}, {5, 5, 5}, 1.0, 0.1, {1.0, 1.0, 1.0});
    // Now check if the cuboid was instantiated as we expect


}


TEST(ForceTest, testForceV1) {
    // Instantiate particles

    // calculate one iteration of the ForceV1 calculation

}


TEST(ForceTest, testForceLennardJones) {
    // Instantiate particles



    // calculate one iteration of the LennardJonesForceIteration

}