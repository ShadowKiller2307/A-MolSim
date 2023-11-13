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

/**
 * A test fixture which initializes some Particles and a container for the following tests
 */
class MolSimTest : public testing::Test {
protected:
    void SetUp() override {
        particles.emplace_back(particleOne);
        particles.emplace_back(particleTwo);
        container.setParticles(particles);
    }

    Particle particleOne{0};
    Particle particleTwo{0};
    std::vector<Particle> particles;
    ParticleContainer container;
    ParticleGenerator particleGenerator;
    ForceV1 forceV1;
    LennardJonesForce lennardJonesForce{5, 1};
};


// unit Test for the particle container
TEST_F(MolSimTest, testGetParticles) {
    EXPECT_EQ(2, container.getParticles()->size());
}

/**
 * Check the values of the particleContainer after the instantiateCuboid method was
 * called
 */
TEST_F(MolSimTest, testGenerateParticlesGenerator) {
    // Instantiate a generator and container for the instantiateCuboid function
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
    //maybe also check the velocities
}

/**
 * Test the ForceV1Calculation against hard coded values
 */
TEST_F(MolSimTest, testForceV1) {
    // calculate one iteration of the ForceV1 calculation
    forceV1.calculateForces(particles);
    // check against hard coded values
    std::array<double, 3> expectedValuesOne{0.0, 0.0, 0.0}; //TODO: durch händisch ausgerechnete Werte ersetzen
    std::array<double, 3> expectedValuesTwo{0.0, 0.0, 0.0}; //TODO: durch händisch ausgerechnete Werte ersetzen
    EXPECT_EQ(particles.at(0).getF(), expectedValuesOne);
    EXPECT_EQ(particles.at(1).getF(), expectedValuesTwo);
}

/**
 * @brief  this test checks whether our old forceV1 calculation and the
 * forceV1 calculation with the lambda generate the same values
 */
TEST_F(MolSimTest, testEqualLambdaForceV1) {
    std::vector<Particle> particlesCopy = particles;
    forceV1.calculateForces(particles);
    forceV1.calculateForcesWithLambda(particlesCopy);
    EXPECT_EQ(particles.at(0).getF(), particlesCopy.at(0).getF());
    EXPECT_EQ(particles.at(1).getF(), particlesCopy.at(1).getF());
}

/**
 * Test the LennardJonesForceCalculation against hard coded values
 */
TEST_F(MolSimTest, testForceLennardJones) {
    // calculate one iteration of the LennardJonesForceIteration
    lennardJonesForce.calculateForces(particles);
    // check against hardcoded values
    std::array<double, 3> expectedValuesOne{0.0, 0.0, 0.0}; //TODO: durch händisch ausgerechnete Werte ersetzen
    std::array<double, 3> expectedValuesTwo{0.0, 0.0, 0.0}; //TODO: durch händisch ausgerechnete Werte ersetzen
    EXPECT_EQ(particles.at(0).getF(), expectedValuesOne);
    EXPECT_EQ(particles.at(1).getF(), expectedValuesTwo);
}

/**
 * @brief  this test checks whether our old lennardjonesforce calculation and the
 * lennardjonesforce calculation with the lambda generate the same values
 */
TEST_F(MolSimTest, testEqualLambdaLennardJonesForce) {
    std::vector<Particle> particlesCopy = particles;
    lennardJonesForce.calculateForces(particles);
    lennardJonesForce.calculateForcesWithLambda(particlesCopy);
    EXPECT_EQ(particles.at(0).getF(), particlesCopy.at(0).getF());
    EXPECT_EQ(particles.at(1).getF(), particlesCopy.at(1).getF());
}