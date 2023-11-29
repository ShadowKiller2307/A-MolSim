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
//#include "../src/xmlSchema/XMLReader.h"
#include "../src/particleContainers/ParticleContainer.h"
#include "../src/particleContainers/ParticleContainerDirSum.h"
#include "../src/particleContainers/ParticleContainerLinCel.h"
#include "../src/boundaryConditions/BoundaryCondition.h"
#include "../src/boundaryConditions/Overflow.h"
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
       containerDirSum.add({2.0, 0.0, 0.0},  {0.0, 0.0, 0.0}, 1, 0);
       containerCuboid = new ParticleContainerDirSum{0.5, 1};
    }
    ParticleContainerDirSum containerDirSum{0.5, 1};
    std::array<double, 3> domainSize{3.0, 3.0, 1.0};
    double cutoffRadius{1.0};
    Reflecting cond1;
    //std::vector<BoundaryCondition> conditions{cond1};
    //ParticleContainerLinCel containerLinCel{0.5, 1, domainSize, cutoffRadius, conditions};   //std::array<double, 3> domainSize, double cutoffRadius, std::vector<BoundaryCondition> &conditions
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
    // TODO: Check the Brownian Motion here, but don't know how yet
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


///This test checks if the one number of cuboids in the xml file is retrieved correctly
/*
TEST_F(MolSimTest, testSimpleCuboid) {
    std::string path = "../../Tests/xmlTestInput/simpleCuboid.xml";

    XMLReader xmlReader(path);

    xmlReader.extractCuboid();
    EXPECT_EQ(1, xmlReader.getCuboidConstructors().size());


}
///This test checks if the parameters of the cuboid are retrieved correctly

TEST_F(MolSimTest, testSimpleCuboidParameters) {
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

    EXPECT_EQ(llfc.at(0),1.3);
    EXPECT_EQ(llfc.at(1),2.0);
    EXPECT_EQ(llfc.at(2),3.0);

    EXPECT_EQ(particlesPerDimension.at(0),15);
    EXPECT_EQ(particlesPerDimension.at(1),20);
    EXPECT_EQ(particlesPerDimension.at(2),30);

    EXPECT_EQ(particleVelocity.at(0),0.4);
    EXPECT_EQ(particleVelocity.at(1),0.4);
    EXPECT_EQ(particleVelocity.at(2),0.4);

    EXPECT_EQ(h,1.89);
    EXPECT_EQ(mass,4.0);
    EXPECT_EQ(type,1);



}
///This test checks if the simulation parameters are extracted correctly
TEST_F(MolSimTest, testSimpleSimulationParameters) {
    std::string path = "../../Tests/xmlTestInput/simpleSimulation.xml";

    XMLReader xmlReader(path);

    xmlReader.extractSimulationParameters();

    SimulationConstructor simulationConstructor = xmlReader.getSimulationConstructor();

    EXPECT_EQ(simulationConstructor.getBaseName(), "DemoSimulation2");
    EXPECT_EQ(simulationConstructor.getWriteFrequency(), 5);
    EXPECT_EQ(simulationConstructor.getT_end(), 101.0);
    EXPECT_EQ(simulationConstructor.getDelta_t(), 0.7);
    EXPECT_EQ(simulationConstructor.getLogLevel(), 2);

    EXPECT_EQ(simulationConstructor.getDomainSize().at(0),120);
    EXPECT_EQ(simulationConstructor.getDomainSize().at(1),50);
    EXPECT_EQ(simulationConstructor.getDomainSize().at(2),1);

    EXPECT_EQ(simulationConstructor.getContainerType(),"LinCel");

}*/
