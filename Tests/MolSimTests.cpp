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
TEST(ParticleContainerTest, testGetParticles)
{
    EXPECT_EQ(2, 2);
    ParticleContainer container;
    std::vector particles{Particle(0), Particle(1)};
    container.setParticles(particles);
    EXPECT_EQ(2, container.getParticles()->size());
}

// Check if the particle generation
/*
TEST(ParticleGeneratorTest, testGenerateParticlesGenerator) {
    // Instantiate a generator and container for the instantiateCuboid function
    ParticleGenerator particleGenerator;
    ParticleContainer container;
    particleGenerator.instantiateCuboid(container, {0.0, 0.0, 0.0}, {2, 2, 2}, 1.0, 0.1, {1.0, 1.0, 1.0}, 0);
    // Now check if the cuboid was instantiated with the particle positions as we expect
    EXPECT_EQ(8, container.getParticles()->size());
    std::array test{0.0, 0.0, 0.0};
    EXPECT_EQ(test, container.getParticles()->at(0).getX());
    test = {0.0, 0.0, 1.0};
    EXPECT_EQ(test, container.getParticles()->at(1).getX());
    test = {0.0, 1.0, 0.0};
    EXPECT_EQ(test, container.getParticles()->at(2).getX());
    test = {0.0, 1.0, 1.0};
    EXPECT_EQ(test, container.getParticles()->at(3).getX());
    test = {1.0, 0.0, 0.0};
    EXPECT_EQ(test, container.getParticles()->at(4).getX());
    test = {1.0, 0.0, 1.0};
    EXPECT_EQ(test, container.getParticles()->at(5).getX());
    test = {1.0, 1.0, 0.0};
    EXPECT_EQ(test, container.getParticles()->at(6).getX());
    test = {1.0, 1.0, 1.0};
    EXPECT_EQ(test, container.getParticles()->at(7).getX());
}


TEST(ForceTest, testForceV1) {
    // Instantiate particles
    // calculate one iteration of the ForceV1 calculation


}


TEST(ForceTest, testForceLennardJones) {
    // Instantiate particles

*/

// calculate one iteration of the LennardJonesForceIteration

//}