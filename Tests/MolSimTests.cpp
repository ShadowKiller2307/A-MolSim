/**
 * @brief Class for Unit Tests which will check the correctness of our code
 */

#include <gtest/gtest.h>
#include <vector>
#include "particleContainers/ParticleContainerDirSum.h"

#include "../src/Particle.h"
#include "../src/particleGenerator/ParticleGenerator.h"
#include "../src/forces/Force.h"
#include "../src/forces/GravPot.h"
#include "../src/forces/LennJon.h"
#include "../src/xmlSchema/XMLReader.h"
#include "../src/particleContainers/ParticleContainer.h"
#include "../src/particleContainers/ParticleContainerDirSum.h"
#include "../src/particleContainers/ParticleContainerLinCel.h"
#include "../src/boundaryConditions/BoundaryCondition.h"
#include "../src/boundaryConditions/Outflow.h"
#include "../src/boundaryConditions/Reflecting.h"

// Difference ASSERT vs EXPECT macros
// ASSERT -> fatal failures
// EXPECT -> nonfatal failures

/**
 * @brief: A test fixture which initializes some Particles and a container for the following tests
 */
class MolSimTest : public testing::Test
{
protected:
    void SetUp() override
    {
        containerDirSum.add({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1, 0);
        containerDirSum.add({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1, 0);
        containerDirSum.add({2.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1, 0);
        linCel2.add({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1, 0);
        linCel2.add({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1, 0);
        linCel2.add({2.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1, 0);
        containerCuboid = new ParticleContainerDirSum{0.5, 1};
    }
    ParticleContainerDirSum containerDirSum{0.5, 1};
    std::array<double, 3> domainSize{3.0, 3.0, 1.0};
    double cutoffRadius{1.0};
    Reflecting *cond1;
    Reflecting *cond2;
    Reflecting *cond3;
    Reflecting *cond4;
    std::vector<BoundaryCondition> conditions{*cond1, *cond2, *cond3, *cond4};
    ParticleContainerLinCel containerLinCel{0.5, 1, domainSize, cutoffRadius, conditions}; // std::array<double, 3> domainSize, double cutoffRadius, std::vector<BoundaryCondition> &conditions
    ParticleContainerLinCel linCel2{0.5, 1, {3.0, 3.0, 1.0}, 1.5, conditions};
    LennJon lennJon{5.0, 1.0};
    GravPot gravPot{};
    ParticleContainer *containerCuboid;
    particleGenerator generator{};
};

/**
 * @brief: very simple test to check whether the container sets the particles correctly
 */
TEST_F(MolSimTest, testGetParticles)
{
    EXPECT_EQ(3, containerDirSum.getParticles().size());
}

/**
 * @brief: Check the position values of the particles in the particleContainer after the instantiateCuboid method was
 * called
 */

TEST_F(MolSimTest, testGenerateParticlesGenerator)
{
    // Instantiate a generator and container for the instantiateCuboid function
    std::array<double, 3> startV{0.0, 0.0, 0.0};
    generator.instantiateCuboid(&containerCuboid, {0.0, 0.0, 0.0}, {2, 2, 2}, startV, 1.0, 1, 0);
    // Now check if the cuboid was instantiated with the particle positions as we expect
    EXPECT_EQ(8, containerCuboid->getParticles().size());
    std::array<double, 3> test{0.0, 0.0, 0.0};
    EXPECT_EQ(test, containerCuboid->getParticles().at(0).getX());
    test = {0.0, 0.0, 1.0};
    EXPECT_EQ(test, containerCuboid->getParticles().at(1).getX());
    test = {0.0, 1.0, 0.0};
    EXPECT_EQ(test, containerCuboid->getParticles().at(2).getX());
    test = {0.0, 1.0, 1.0};
    EXPECT_EQ(test, containerCuboid->getParticles().at(3).getX());
    test = {1.0, 0.0, 0.0};
    EXPECT_EQ(test, containerCuboid->getParticles().at(4).getX());
    test = {1.0, 0.0, 1.0};
    EXPECT_EQ(test, containerCuboid->getParticles().at(5).getX());
    test = {1.0, 1.0, 0.0};
    EXPECT_EQ(test, containerCuboid->getParticles().at(6).getX());
    test = {1.0, 1.0, 1.0};
    EXPECT_EQ(test, containerCuboid->getParticles().at(7).getX());
    //  Now check if the sphere was instantiated with the particle positions as we expect
    test = {4.0, 4.0, 0.0};
    EXPECT_EQ(test, containerCuboid->getParticles().at(8).getX());
    test = {4.0, 5.0, 0.0};
    EXPECT_EQ(test, containerCuboid->getParticles().at(9).getX());
    test = {4.0, 6.0, 0.0};
    EXPECT_EQ(test, containerCuboid->getParticles().at(10).getX());
    test = {5.0, 4.0, 0.0};
    EXPECT_EQ(test, containerCuboid->getParticles().at(11).getX());
    test = {5.0, 5.0, 0.0};
    EXPECT_EQ(test, containerCuboid->getParticles().at(12).getX());
    test = {5.0, 6.0, 0.0};
    EXPECT_EQ(test, containerCuboid->getParticles().at(13).getX());
    test = {6.0, 4.0, 0.0};
    EXPECT_EQ(test, containerCuboid->getParticles().at(14).getX());
    test = {6.0, 5.0, 0.0};
    EXPECT_EQ(test, containerCuboid->getParticles().at(15).getX());
    test = {6.0, 6.0, 0.0};
    EXPECT_EQ(test, containerCuboid->getParticles().at(16).getX());
    // Check for two random particles(one of the sphere and one of the cuboid) if the mass was set correctly
    EXPECT_EQ(1.0, containerCuboid->getParticles().at(4).getM());
    EXPECT_EQ(1.0, containerCuboid->getParticles().at(16).getM());
}

/**
 * @brief: Generate a cuboid and a sphere for the linked cells
 */
TEST_F(MolSimTest, testGenerateParticlesLinCelContainer)
{
    ParticleContainer *cuboidLinkedCel = &containerLinCel;
    std::array<double, 3> startV{0.0, 0.0, 0.0};
    generator.instantiateCuboid(&cuboidLinkedCel, {0.5, 0.5, 0.0}, {2, 2, 0}, startV, 1.0, 1, 0);
    for (int i = 0; i < containerLinCel.getAmountOfCells(); ++i)
    {
        std::cout << containerLinCel.getCells().at(i).size() << std::endl;
    }

    /*EXPECT_EQ(containerLinCel.getAmountOfCells(), 25);
    std::array<double, 3> test{0.5, 0.5, 0.0};
    EXPECT_EQ(test, containerLinCel.getCells().at(6).at(0).getX());
    test = {0.5, 1.5, 0.0};
    EXPECT_EQ(test, containerLinCel.getCells().at(11).at(0).getX());
    test = {1.5, 0.5, 0.0};
    EXPECT_EQ(test, containerLinCel.getCells().at(7).at(0).getX());
    test = {1.5, 1.5, 0.0};
    EXPECT_EQ(test, containerLinCel.getCells().at(12).at(0).getX());*/
}

/**
 * @brief: check force calculation for Lennard Jones for the Linked cells
 */
TEST_F(MolSimTest, testForcesLinkedCells)
{
    linCel2.setForce(lennJon.innerPairs());
    std::cout << "Calculate the forces" << std::endl;
    linCel2.calculateForces();
    EXPECT_EQ(linCel2.getAmountOfCells(), 16);
    EXPECT_EQ(linCel2.getAmountOfParticles(), 3);
    for (int i = 0; i < linCel2.getAmountOfCells(); ++i)
    {
        if (linCel2.getCells().at(i).size() > 0)
        {
            std::cout << i << linCel2.getCells().at(i).at(0) << std::endl;
            if (linCel2.getCells().at(i).size() > 1)
            {
                std::cout << i << linCel2.getCells().at(i).at(1) << std::endl;
            }
        }
    }
    // check against hardcoded values
    std::array<double, 3> expectedValuesOne{-119.091796875, 0.0, 0.0};
    std::array<double, 3> expectedValuesTwo{0.0, 0.0, 0.0};
    std::array<double, 3> expectedValuesThree{119.091796875, 0.0, 0.0};
    double test = (465.0 / 256.0);
    std::cout << test << std::endl;
    // printf("Test=%.17le", test);
    EXPECT_EQ(containerDirSum.getParticles().at(0).getF(), expectedValuesOne);
    EXPECT_EQ(containerDirSum.getParticles().at(1).getF(), expectedValuesTwo);
    EXPECT_EQ(containerDirSum.getParticles().at(2).getF(), expectedValuesThree);
}

/**
 * @brief: check if a single Particle in a Boundary Cell, that moves towards the border of the domain,
 * stays within the domain when the Boundary is set to Reflecting
 */
TEST_F(MolSimTest, testReflectingBoundary)
{
    // left domain border should have the BoundaryCondition Reflecting
    containerLinCel.add({0.5, 1.5, 0.0}, {-1.0, 0.0, 0.0}, 1, 0);
    for (int i = 0; i < 10; ++i)
    {
        // check whether the particle leaves the domain and gets deleted
        //  calculate new x
        containerLinCel.calculatePosition();
        // calculate new f
        containerLinCel.iterBoundary();
        containerLinCel.calculateForces();
        // calculate new v
        containerLinCel.calculateVelocity();
        EXPECT_EQ(containerLinCel.getAmountOfParticles(), 1);
    }
}

/**
 *  @brief: check if a single Particle in a Boundary Cell, that moves towards the border of the domain,
 *  gets deleted when leaving the cell if the Boundary is set to Overflow
 */
TEST_F(MolSimTest, testOverflowBoundary)
{
    // left domain border should have the BoundaryCondition Overflow
    containerLinCel.add({0.5, 1.5, 0.0}, {-1.0, 0.0, 0.0}, 1, 0);
    for (int i = 0; i < 10; ++i)
    {
        // check whether the particle leaves the domain and gets deleted
        //  calculate new x
        containerLinCel.calculatePosition();
        // calculate new f
        containerLinCel.iterHalo();
        containerLinCel.calculateForces();
        // calculate new v
        containerLinCel.calculateVelocity();
    }
    EXPECT_EQ(containerLinCel.getAmountOfParticles(), 0);
}

/**
 * @brief: Test the ForceV1Calculation against hard coded values
 */

TEST_F(MolSimTest, testForceV1)
{
    containerDirSum.setForce(gravPot.innerPairs());
    containerDirSum.calculateForces();
    // check against hard coded values
    std::array<double, 3> expectedValuesOne{1.25, 0.0, 0.0};
    std::array<double, 3> expectedValuesTwo{0.0, 0.0, 0.0};
    std::array<double, 3> expectedValuesThree{-1.25, 0.0, 0.0};
    EXPECT_EQ(containerDirSum.getParticles().at(0).getF(), expectedValuesOne);
    EXPECT_EQ(containerDirSum.getParticles().at(1).getF(), expectedValuesTwo);
    EXPECT_EQ(containerDirSum.getParticles().at(2).getF(), expectedValuesThree);
}

/**
 * @brief: Test the LennardJonesForceCalculation against hard coded values
 */

TEST_F(MolSimTest, testForceLennardJones)
{
    // calculate one iteration of the LennardJonesForceIteration
    containerDirSum.setForce(lennJon.innerPairs());
    containerDirSum.calculateForces();
    // check against hardcoded values
    std::array<double, 3> expectedValuesOne{-119.091796875, 0.0, 0.0};
    // have to check whether this is due to the double, an error in the calculation in the program or an
    // error in the calculation on paper
    std::array<double, 3> expectedValuesTwo{0.0, 0.0, 0.0};
    std::array<double, 3> expectedValuesThree{119.091796875, 0.0, 0.0};
    double test = (465.0 / 256.0);
    std::cout << test << std::endl;
    printf("Test=%.17le", test);
    EXPECT_EQ(containerDirSum.getParticles().at(0).getF(), expectedValuesOne);
    EXPECT_EQ(containerDirSum.getParticles().at(1).getF(), expectedValuesTwo);
    EXPECT_EQ(containerDirSum.getParticles().at(2).getF(), expectedValuesThree);
}

/// This test checks if the one number of cuboids in the xml file is retrieved correctly

TEST_F(MolSimTest, testSimpleCuboid)
{
    std::string path = "../../Tests/xmlTestInput/simpleCuboid.xml";

    XMLReader xmlReader(path);

    xmlReader.extractCuboid();
    EXPECT_EQ(1, xmlReader.getCuboidConstructors().size());
}
/// This test checks if the parameters of the cuboid are retrieved correctly

TEST_F(MolSimTest, testSimpleCuboidParameters)
{
    std::string path = "../../Tests/xmlTestInput/simpleCuboid.xml";

    XMLReader xmlReader(path);

    xmlReader.extractCuboid();
    auto cuboidConstructor = xmlReader.getCuboidConstructors().at(0);
    auto llfc = cuboidConstructor.getLlfc();
    auto particlesPerDimension = cuboidConstructor.getParticlesPerDimension();
    auto particleVelocity = cuboidConstructor.getParticleVelocity();
    double h = cuboidConstructor.getH();
    double mass = cuboidConstructor.getMass();
    int type = cuboidConstructor.getType();

    EXPECT_EQ(llfc.at(0), 1.3);
    EXPECT_EQ(llfc.at(1), 2.0);
    EXPECT_EQ(llfc.at(2), 3.0);

    EXPECT_EQ(particlesPerDimension.at(0), 15);
    EXPECT_EQ(particlesPerDimension.at(1), 20);
    EXPECT_EQ(particlesPerDimension.at(2), 30);

    EXPECT_EQ(particleVelocity.at(0), 0.4);
    EXPECT_EQ(particleVelocity.at(1), 0.4);
    EXPECT_EQ(particleVelocity.at(2), 0.4);

    EXPECT_EQ(h, 1.89);
    EXPECT_EQ(mass, 4.0);
    EXPECT_EQ(type, 1);
}
/// This test checks if the simulation parameters are extracted correctly
TEST_F(MolSimTest, testSimpleSimulationParameters)
{
    std::string path = "../../Tests/xmlTestInput/simpleSimulation.xml";

    XMLReader xmlReader(path);

    xmlReader.extractSimulationParameters();

    SimulationConstructor simulationConstructor = xmlReader.getSimulationConstructor();

    EXPECT_EQ(simulationConstructor.getBaseName(), "DemoSimulation2");
    EXPECT_EQ(simulationConstructor.getWriteFrequency(), 5);
    EXPECT_EQ(simulationConstructor.getT_end(), 101.0);
    EXPECT_EQ(simulationConstructor.getDelta_t(), 0.7);
    EXPECT_EQ(simulationConstructor.getLogLevel(), 2);

    EXPECT_EQ(simulationConstructor.getDomainSize().at(0), 120);
    EXPECT_EQ(simulationConstructor.getDomainSize().at(1), 50);
    EXPECT_EQ(simulationConstructor.getDomainSize().at(2), 1);

    EXPECT_EQ(simulationConstructor.getContainerType(), "LinCel");
}
