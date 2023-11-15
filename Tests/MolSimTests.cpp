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
        particles.emplace_back(particleThree);
        container.setParticles(particles);
    }

    std::array<double,3> temp{0.0, 0.0, 0.0};
    std::array<double,3> temp2{1.0, 0.0, 0.0};
    std::array<double, 3> temp3{2.0, 0.0, 0.0};
    std::array<double,3> temp4{0.0, 0.0, 0.0};
    std::array<double,3> temp5{1.0, 0.0, 0.0};
    std::array<double, 3> temp6{2.0, 0.0, 0.0};
    Particle particleOne{temp, temp, 1, 0};
    Particle particleTwo{temp2, temp, 1, 0};
    Particle particleThree{temp3, temp, 1, 0};
    std::vector<Particle> particles;
    ParticleContainer container;
    ParticleContainer container2;
    ParticleGenerator particleGenerator;
    ForceV1 forceV1;
    LennardJonesForce lennardJonesForce{5, 1};
};


// unit Test for the particle container
TEST_F(MolSimTest, testGetParticles) {
    EXPECT_EQ(2, container.getParticles().size());
}

/**
 * Check the values of the particleContainer after the instantiateCuboid method was
 * called
 */
TEST_F(MolSimTest, testGenerateParticlesGenerator) {
    // Instantiate a generator and container for the instantiateCuboid function
    std::array<double, 3> startV{0.0, 0.0, 0.0};
    ParticleContainer containerCuboid;
    particleGenerator.instantiateCuboid(containerCuboid, {0.0, 0.0, 0.0}, {2, 2, 2}, startV, 1.0, 1, 0);
    // Now check if the cuboid was instantiated with the particle positions as we expect
    EXPECT_EQ(8, containerCuboid.getParticles().size());
    std::array test{0.0, 0.0, 0.0};
    EXPECT_EQ(test, containerCuboid.getParticles().at(0).getX());
    test = {0.0, 0.0, 1.0};
    EXPECT_EQ(test, containerCuboid.getParticles().at(1).getX());
    test = {0.0, 1.0, 0.0};
    EXPECT_EQ(test, containerCuboid.getParticles().at(2).getX());
    test = {0.0, 1.0, 1.0};
    EXPECT_EQ(test, containerCuboid.getParticles().at(3).getX());
    test = {1.0, 0.0, 0.0};
    EXPECT_EQ(test, containerCuboid.getParticles().at(4).getX());
    test = {1.0, 0.0, 1.0};
    EXPECT_EQ(test, containerCuboid.getParticles().at(5).getX());
    test = {1.0, 1.0, 0.0};
    EXPECT_EQ(test, containerCuboid.getParticles().at(6).getX());
    test = {1.0, 1.0, 1.0};
    EXPECT_EQ(test, containerCuboid.getParticles().at(7).getX());
    //maybe also check the velocities
}

/**
 * Test the ForceV1Calculation against hard coded values
 */
TEST_F(MolSimTest, testForceV1) {
    // calculate one iteration of the ForceV1 calculation
    forceV1.calculateForces(particles);
    // check against hard coded values
    std::array<double, 3> expectedValuesOne{1.25, 0.0, 0.0};
    std::array<double, 3> expectedValuesTwo{0.0, 0.0, 0.0};
    std::array<double, 3> expectedValuesThree{-1.25, 0.0, 0.0};
    EXPECT_EQ(particles.at(0).getF(), expectedValuesOne);
    EXPECT_EQ(particles.at(1).getF(), expectedValuesTwo);
    EXPECT_EQ(particles.at(2).getF(), expectedValuesThree);
}

/**
 * @brief  this test checks whether our old forceV1 calculation and the
 * forceV1 calculation with the lambda generate the same values
 */
TEST_F(MolSimTest, testEqualLambdaForceV1) {
    auto particlesCopy = particles;
    container2.setParticles(particlesCopy);
    forceV1.calculateForces(particles);
    forceV1.calculateForcesWithLambda(container2);
    EXPECT_EQ(particles.at(0).getF(), container2.getParticles().at(0).getF());
    EXPECT_EQ(particles.at(1).getF(), container2.getParticles().at(1).getF());
    EXPECT_EQ(particles.at(2).getF(), container2.getParticles().at(2).getF());
}

/**
 * Test the LennardJonesForceCalculation against hard coded values
 */
TEST_F(MolSimTest, testForceLennardJones) {
    // calculate one iteration of the LennardJonesForceIteration
    lennardJonesForce.calculateForces(particles);
    // check against hardcoded values
    std::array<double, 3> expectedValuesOne{-119.091796875, 0.0, 0.0};
    //have to check whether this is due to the double, an error in the calculation in the program or an
    //error in the calculation on paper
    std::array<double, 3> expectedValuesTwo{0.0, 0.0, 0.0};
    std::array<double, 3> expectedValuesThree{119.091796875, 0.0, 0.0};
    double test = (465.0/256.0);
    std::cout << test << std::endl;
    printf("Test=%.17le", test);
    EXPECT_EQ(particles.at(0).getF(), expectedValuesOne);
    EXPECT_EQ(particles.at(1).getF(), expectedValuesTwo);
    EXPECT_EQ(particles.at(2).getF(), expectedValuesThree);
    }

/**
 * @brief  this test checks whether our old lennardjonesforce calculation and the
 * lennardjonesforce calculation with the lambda generate the same values
 */
TEST_F(MolSimTest, testEqualLambdaLennardJonesForce) {
    auto particlesCopy = particles;
    container2.setParticles(particlesCopy);
    lennardJonesForce.calculateForces(particles);
    lennardJonesForce.calculateForcesWithLambda(container2);
    EXPECT_EQ(particles.at(0).getF(), container2.getParticles().at(0).getF());
    EXPECT_EQ(particles.at(1).getF(), container2.getParticles().at(1).getF());
    EXPECT_EQ(particles.at(2).getF(), container2.getParticles().at(2).getF());
}
