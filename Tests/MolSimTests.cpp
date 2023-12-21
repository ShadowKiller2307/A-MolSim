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
        containerDirSum2.add({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1, 0);
        containerDirSum2.add({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1, 0);
        containerDirSum2.add({2.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1, 0);
        linCel2.add({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1, 0);
        linCel2.add({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1, 0);
        linCel2.add({2.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1, 0);
        containerCuboid = new ParticleContainerDirSum{0.5, 1, 1};
    }
    std::array<double, 3> domainSize{3.0, 3.0, 1.0};
    double cutoffRadius{1.0};
    LennJon lennJon{5.0, 1.0};
    GravPot gravPot1{};
    ParticleContainerDirSum containerDirSum{0.5, 1, 1};
    ParticleContainerDirSum containerDirSum2{0.5, 1, 1};
    ParticleContainerLinCel containerLinCel{0.5, 1, 1, domainSize, "rrrrrr", 1.0}; // std::array<double, 3> domainSize, double cutoffRadius, std::vector<BoundaryCondition> &conditions
    ParticleContainerLinCel linCel2{0.5, 1, 1, {3.0, 3.0, 1.0}, "oooooo", 1.5};
    ParticleContainerLinCel linCel3{0.5, 1, 1, {5.0, 5.0, 5.0}, "oooooo", 1.0};
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
/*
ParticleContainerLinCel linCelInitTemp{0.001, 10, 1, {120.0, 120.0, 1.0}, "rrrrrr", 3.0,
                                       true, 10, true, 40, 40, 1.0, 0}; // parameters for the thermostat
*/

/*
void particleGenerator::instantiateCuboid(ParticleContainer **container, const std::array<double, 3> &llfc,
                                          const std::array<unsigned int, 3> &particlesPerDimension,
                                          std::array<double, 3> &particleVelocity, double h, double m, int type, double epsilon, double sigma, double initT = 0.0)
*/


/*void particleGenerator::instantiateSphere(ParticleContainer **container, const std::array<double, 3> &center,
                                          const int32_t &sphereRadius, const std::array<double, 3> &particleVelocity,
                                          double h, double m, bool is2D, int type, double epsilon, double sigma, double initT = 0.0)
{*/
TEST_F(MolSimTest, testGenerateParticlesGenerator)
{
    // Instantiate a generator and container for the instantiateCuboid function
    std::array<double, 3> startV{0.0, 0.0, 0.0};
    particleGenerator::instantiateCuboid(&containerCuboid, {0.0, 0.0, 0.0}, {2, 2, 2}, startV, 1.0, 1, 0, 5, 1, 40);
    particleGenerator::instantiateSphere(&containerCuboid, {5.0, 5.0, 0.0}, 2, startV, 1.0, 1, true, 1, 5, 1, 40);
    // Now check if the cuboid was instantiated with the particle positions as we expect
    EXPECT_EQ(17, containerCuboid->getParticles().size());
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


TEST_F(MolSimTest, testPositionIndexTranslation)
{
    std::array<double, 3> pos = {0.5, 0.5, 0};
    EXPECT_EQ(31, containerLinCel.translate3DPosTo1D(pos));
    std::array<double, 3> pos2 = {-0.5, -0.5, 0};
    EXPECT_EQ(25, containerLinCel.translate3DPosTo1D(pos2));
}

TEST_F(MolSimTest, test3DIndexTo1DIndexTranslation)
{
    int index = containerLinCel.translate3DIndTo1D(0, 0, 0);
    EXPECT_EQ(0, index);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief: Generate a cuboid and a sphere for the linked cells
 */

/*
void particleGenerator::instantiateCuboid(ParticleContainer **container, const std::array<double, 3> &llfc,
                                          const std::array<unsigned int, 3> &particlesPerDimension,
                                          std::array<double, 3> &particleVelocity, double h, double m, int type, double epsilon, double sigma, double initT = 0.0)
*/
TEST_F(MolSimTest, testGenerateParticlesLinCelContainer)
{
    ParticleContainer *cuboidLinkedCel = &containerLinCel;
    std::array<double, 3> startV{0.0, 0.0, 0.0};
    particleGenerator::instantiateCuboid(&cuboidLinkedCel, {0.5, 0.5, 0.0}, {2, 2, 1}, startV, 1.0, 1, 0, 5, 1, 40);
    std::cout << containerLinCel.getCells().size() << std::endl;
    EXPECT_EQ(containerLinCel.getAmountOfCells(), 75);
    std::array<double, 3> test{0.5, 0.5, 0.0};
    EXPECT_EQ(test, containerLinCel.getCells().at(31).at(0).getX());
    test = {0.5, 1.5, 0.0};
    EXPECT_EQ(test, containerLinCel.getCells().at(36).at(0).getX());
    test = {1.5, 0.5, 0.0};
    EXPECT_EQ(test, containerLinCel.getCells().at(32).at(0).getX());
    test = {1.5, 1.5, 0.0};
    EXPECT_EQ(test, containerLinCel.getCells().at(37).at(0).getX());
}

/**
 * @brief: check if a single Particle in a Boundary Cell, that moves towards the border of the domain,
 * stays within the domain when the Boundary is set to Reflecting
 */
TEST_F(MolSimTest, testReflectingBoundary)
{
    // left domain border should have the BoundaryCondition Reflecting
    containerLinCel.add({0.5, 1.5, 0.0}, {-1.0, 0.0, 0.0}, 1, 0);
    /*EXPECT_EQ(containerLinCel.affectsForce(0), true);
    EXPECT_EQ(containerLinCel.affectsForce(5), false);*/
    for (int i = 0; i < 10; ++i)
    {
        // check whether the particle leaves the domain and gets deleted
        //  calculate new x
        containerLinCel.calculatePosition();
        // calculate new f
        containerLinCel.calculateForces();
        // calculate new v
        containerLinCel.calculateVelocity();
        EXPECT_EQ(containerLinCel.getAmountOfParticles(), 1);
    }
}

/**
 * @brief: check force calculation for Lennard Jones for the Linked cells
 */
TEST_F(MolSimTest, testForcesLinkedCells)
{
    std::cout << "Calculate the forces" << std::endl;
    linCel2.calculateForces();
    //  EXPECT_EQ(linCel2.getAmountOfCells(), 16);
    //  EXPECT_EQ(linCel2.getAmountOfParticles(), 3);
    // check against hardcoded values
    // different values than the direct sum container as the distance between the particle on the left and the particle on the right
    // are bigger than the cutoff radius
    std::array<double, 3> expectedValuesOne{-120, 0.0, 0.0};
    std::array<double, 3> expectedValuesTwo{0.0, 0.0, 0.0};
    std::array<double, 3> expectedValuesThree{120, 0.0, 0.0};
    EXPECT_EQ(linCel2.getCells().at(21).at(0).getF(), expectedValuesOne);
    EXPECT_EQ(linCel2.getCells().at(21).at(1).getF(), expectedValuesTwo);
    EXPECT_EQ(linCel2.getCells().at(22).at(0).getF(), expectedValuesThree);
}

/**
 *  @brief: check if a single Particle in a Boundary Cell, that moves towards the border of the domain,
 *  gets deleted when leaving the cell if the Boundary is set to Overflow
 */
/*TEST_F(MolSimTest, testOverflowBoundary)
{
    // left domain border should have the BoundaryCondition Overflow
    ParticleContainerLinCel containerLinCel2{0.5, 1, 1, domainSize, "oooooo", lennJon, 1.0};
    // containerLinCel.add({0.5, 1.5, 0.0}, {-1.0, 0.0, 0.0}, 1, 0);
    for (int i = 0; i < 10; ++i)
    {
        // check whether the particle leaves the domain and gets deleted
        //  calculate new x
        containerLinCel2.calculatePosition();
        // calculate new f
        containerLinCel2.calculateForces();
        // calculate new v
        containerLinCel2.calculateVelocity();
    }
    EXPECT_EQ(containerLinCel2.getAmountOfParticles(), 0);
}*/

/**
 * @brief: Test the ForceV1Calculation against hard coded values
 */
/*TEST_F(MolSimTest, testForceV1)
{
    containerDirSum2.calculateForces();
    // check against hard coded values
    std::array<double, 3> expectedValuesOne{1.25, 0.0, 0.0};
    std::array<double, 3> expectedValuesTwo{0.0, 0.0, 0.0};
    std::array<double, 3> expectedValuesThree{-1.25, 0.0, 0.0};
    EXPECT_EQ(containerDirSum2.getParticles().at(0)->getF(), expectedValuesOne);
    EXPECT_EQ(containerDirSum2.getParticles().at(1)->getF(), expectedValuesTwo);
    EXPECT_EQ(containerDirSum2.getParticles().at(2)->getF(), expectedValuesThree);
}*/

/**
 * @brief: Test the LennardJonesForceCalculation against hard coded values
 */
/*TEST_F(MolSimTest, testForceLennardJones)
{
    // calculate one iteration of the LennardJonesForceIteration
    containerDirSum.calculateForces();
    // check against hardcoded values
    std::array<double, 3> expectedValuesOne{-119.091796875, 0.0, 0.0};
    // have to check whether this is due to the double, an error in the calculation in the program or an
    // error in the calculation on paper
    std::array<double, 3> expectedValuesTwo{0.0, 0.0, 0.0};
    std::array<double, 3> expectedValuesThree{119.091796875, 0.0, 0.0};
    *//*double test = (465.0 / 256.0);
    std::cout << test << std::endl;
    printf("Test=%.17le", test);*//*
    EXPECT_EQ(containerDirSum.getParticles().at(0)->getF(), expectedValuesOne);
    EXPECT_EQ(containerDirSum.getParticles().at(1)->getF(), expectedValuesTwo);
    EXPECT_EQ(containerDirSum.getParticles().at(2)->getF(), expectedValuesThree);
}*/

////////////////////////////////////////////////////// TESTS FOR THE XML INPUT BEGIN ///////////////////////////////////

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
    double sigma = cuboidConstructor.getSigma();
    double epsilon = cuboidConstructor.getEpsilon();

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
    EXPECT_EQ(sigma,1);
    EXPECT_EQ(epsilon,5);
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
    EXPECT_EQ(simulationConstructor.getBoundaries(), "rororr");
    EXPECT_EQ(simulationConstructor.getCutOffRadius(), 3);
    EXPECT_EQ(simulationConstructor.getGrav(),10);
    EXPECT_EQ(simulationConstructor.isUseThermostat(),true);
    EXPECT_EQ(simulationConstructor.getIsGradual(),true);

    EXPECT_EQ(simulationConstructor.getInitialTemp(),40);
    EXPECT_EQ(simulationConstructor.getNThermostat(),10);
    EXPECT_EQ(simulationConstructor.getTempTarget(),60);
    EXPECT_EQ(simulationConstructor.getMaxDiff(),15);
}
/// This test checks if the parameters of two spheres are retrieved correctly
TEST_F(MolSimTest, testSimpleSphereParameters)
{
    std::string path = "../../Tests/xmlTestInput/simpleSphere.xml";

    XMLReader xmlReader(path);

    xmlReader.extractSphere();

    EXPECT_EQ(xmlReader.getSphereConstructors().size(), 2);

    SphereConstructor sphereConstructor = xmlReader.getSphereConstructors().at(0);
    SphereConstructor sphereConstructor1 = xmlReader.getSphereConstructors().at(1);

    EXPECT_EQ(sphereConstructor.getCenterCoordinates().at(0), 60);
    EXPECT_EQ(sphereConstructor.getCenterCoordinates().at(1), 25);
    EXPECT_EQ(sphereConstructor.getCenterCoordinates().at(2), 0);

    EXPECT_EQ(sphereConstructor.getInitialVelocity().at(0), 0);
    EXPECT_EQ(sphereConstructor.getInitialVelocity().at(1), -10);
    EXPECT_EQ(sphereConstructor.getInitialVelocity().at(2), 0);

    EXPECT_EQ(sphereConstructor.getRadius(), 15);
    EXPECT_EQ(sphereConstructor.getDistance(), 1.1225);
    EXPECT_EQ(sphereConstructor.getMass(), 1.0);
    EXPECT_EQ(sphereConstructor.getType(),0);

    EXPECT_EQ(sphereConstructor1.getCenterCoordinates().at(0), 90);
    EXPECT_EQ(sphereConstructor1.getCenterCoordinates().at(1), 20);
    EXPECT_EQ(sphereConstructor1.getCenterCoordinates().at(2), 0);

    EXPECT_EQ(sphereConstructor1.getInitialVelocity().at(0), 0);
    EXPECT_EQ(sphereConstructor1.getInitialVelocity().at(1), 20);
    EXPECT_EQ(sphereConstructor1.getInitialVelocity().at(2), 0);

    EXPECT_EQ(sphereConstructor1.getRadius(), 19);
    EXPECT_EQ(sphereConstructor1.getDistance(), 1.1226);
    EXPECT_EQ(sphereConstructor1.getMass(), 1.0);
    EXPECT_EQ(sphereConstructor1.getType(),1);
}

//////////////////////////////////////////////// TESTS FOR THE XML INPUT END ///////////////////////////////////////////


/////////////////////////////////////////////// TEST FOR THE BOUNDARY CONDITIONS BEGIN /////////////////////////////////

TEST_F(MolSimTest, testCreateReflecting)
{
    ParticleContainerLinCel linCelTest{0.5, 1, 1, {3.0, 3.0, 1.0}, "rrrrrr", 1.5};
    linCelTest.add({0.5, 1.5, 0.0}, {0.0, 0.0, 0.0}, 1, 0);
    auto lambda = linCelTest.createReflectingLambdaBoundary(0, 0);
    // auto lambda2 = linCelTest.createReflectingLambdaBoundary(1, 0);
    lambda(1, 2, 1);
    // lambda2(1, 2, 1);
    EXPECT_EQ(1, 1);
}

// we used the following tests to debug our reflecting boundaries
/*
TEST_F(MolSimTest, testReflectingOneParticleHorizontalMovement)
{
    ParticleContainerLinCel linCelTest{0.1, 50, 1, {3.0, 3.0, 1.0}, "rrrrrr", 1.0};
    linCelTest.add({1.5, 1.5, 0.5}, {-0.1, 0.0, 0.0}, 1, 0);
    std::unique_ptr<std::array<double, 3>> position = std::make_unique<std::array<double, 3>>(linCelTest.getParticles().at(0).getX());
    std::unique_ptr<std::array<double, 3>> force = std::make_unique<std::array<double, 3>>(linCelTest.getParticles().at(0).getF());
    std::cout << "Index: " << linCelTest.getCells().at(36).size() << std::endl;
    // linCelTest.getCells().at()
    linCelTest.simulateParticles();
}
*/

/*
TEST_F(MolSimTest, testReflectingOneParticleHorizontalMovementRightBoundary)
{
    ParticleContainerLinCel linCelTest{0.1, 50, 1, {3.0, 3.0, 1.0}, "rrrrrr", 1.0};
    linCelTest.add({1.5, 1.5, 0.5}, {0.1, 0.0, 0.0}, 1, 0);
    std::unique_ptr<std::array<double, 3>> position = std::make_unique<std::array<double, 3>>(linCelTest.getParticles().at(0).getX());
    std::unique_ptr<std::array<double, 3>> force = std::make_unique<std::array<double, 3>>(linCelTest.getParticles().at(0).getF());
    std::cout << "Index: " << linCelTest.getCells().at(36).size() << std::endl;
    // linCelTest.getCells().at()
    linCelTest.simulateParticles();
}

TEST_F(MolSimTest, testReflectingOneParticleVerticalMovement)
{
    ParticleContainerLinCel linCelTest{0.01, 50, 1, {3.0, 3.0, 1.0}, "rrrrrr", 1.0};
    linCelTest.add({1.5, 1.5, 0.5}, {0.0, -0.1, 0.0}, 1, 0);
    std::unique_ptr<std::array<double, 3>> position = std::make_unique<std::array<double, 3>>(linCelTest.getParticles().at(0).getX());
    std::unique_ptr<std::array<double, 3>> force = std::make_unique<std::array<double, 3>>(linCelTest.getParticles().at(0).getF());
    // linCelTest.getCells().at()
    linCelTest.simulateParticles();
    std::cout << "Particles: " << linCelTest.getParticles().at(0) << std::endl;
    std::cout << "Temperatur: " << linCelTest.calculateTemperature() << std::endl;
}*/

///////////////////////////////////////////////////// TESTS FOR THE BOUNDARY CONDITIONS END ////////////////////////////


///////////////////////////////////////////////////// TESTS FOR THE THERMOSTAT BEGIN ///////////////////////////////////

/**
 * @brief: In this test we check whether initializing a system with an initial temperature works
 */
TEST_F(MolSimTest, initializeTemperature) {
    // first call the particle generator with the desired initial temperature
    double initialTemp = 40;
    ParticleContainerLinCel linCelInitTemp{0.01, 50, 1, {120.0, 120.0, 1.0}, "rrrrrr", 3.0};
    ParticleContainer *ptr = &linCelInitTemp;
    ParticleContainer **ptrptr = &ptr;
    std::array<double, 3> vel = {0.0, 0.0, 0.0};
    particleGenerator::instantiateCuboid(ptrptr, {1.0, 1.0, 0.0}, {10, 10, 1}, vel, 1, 1.0, 0, 5, 1, 40);
   // double kineticEnergy = linCelInitTemp.calculateKinEnergy();
    double initTemp = linCelInitTemp.calculateTemperature();
    std::cout << initTemp << std::endl;
    //EXPECT_DOUBLE_EQ(initTemp,initialTemp);
    EXPECT_NEAR(initTemp, initialTemp, 5.0);
}

/**
 * @brief This tests checks whether after initializing the system with a temperature around 40 degrees,
 * the heating process of the system works with the thermostat
 * Therefore it should be heated up to around 60 degrees
 */
TEST_F(MolSimTest, testHeating) {
    double initialTemp = 40;
    double targetTemp = 60;
    ParticleContainerLinCel linCelInitTemp{0.001, 10, 1, {120.0, 120.0, 1.0}, "rrrrrr", 3.0,
                                           true, 10, true, 40, 60, 1.0, 0}; // parameters for the thermostat
    ParticleContainer *ptr = &linCelInitTemp;
    ParticleContainer **ptrptr = &ptr;
    std::array<double, 3> vel = {0.0, 0.0, 0.0};
    particleGenerator::instantiateCuboid(ptrptr, {1.0, 1.0, 0.0},
                                         {10, 10, 1}, vel, 1, 1.0, 0, 5, 1, initialTemp);
    EXPECT_NEAR(linCelInitTemp.calculateTemperature(), initialTemp, 5.0);
    linCelInitTemp.simulateParticles();
    EXPECT_NEAR(linCelInitTemp.calculateTemperature(), targetTemp, 5.0);
}

/**
 * @brief This tests checks whether after initializing the systen with a temperature around 40 degrees,
 * the cooling process of the system works with the thermostat
 * Therefore it should be cooled down to around 20 degrees
 */
TEST_F(MolSimTest, testCooling) {
    double initialTemp = 40;
    double targetTemp = 20;
    ParticleContainerLinCel linCelInitTemp{0.001, 10, 1, {120.0, 120.0, 1.0}, "rrrrrr", 3.0,
                                           true, 10, true, 40, 20, 1.0, 0}; // parameters for the thermostat
    ParticleContainer *ptr = &linCelInitTemp;
    ParticleContainer **ptrptr = &ptr;
    std::array<double, 3> vel = {0.0, 0.0, 0.0};
    particleGenerator::instantiateCuboid(ptrptr, {1.0, 1.0, 0.0},
                                         {10, 10, 1}, vel, 1, 1.0, 0, 5, 1, initialTemp);
    EXPECT_NEAR(linCelInitTemp.calculateTemperature(), initialTemp, 5.0);
    linCelInitTemp.simulateParticles();
    EXPECT_NEAR(linCelInitTemp.calculateTemperature(), targetTemp, 5.0);
}

/**
 * @brief: This test checks if the system keeps the initial temperature with a max error of 5.0
 */
TEST_F(MolSimTest, testHoldingATemperature) {
    double initialTemp = 40;
    double targetTemp = 40;
    ParticleContainerLinCel linCelInitTemp{0.001, 10, 1, {120.0, 120.0, 1.0}, "rrrrrr", 3.0,
                                           true, 10, true, 40, 40, 1.0, 0}; // parameters for the thermostat
    ParticleContainer *ptr = &linCelInitTemp;
    ParticleContainer **ptrptr = &ptr;
    std::array<double, 3> vel = {0.0, 0.0, 0.0};
    // initialize a cuboid with an initial temperature of 40
    particleGenerator::instantiateCuboid(ptrptr, {1.0, 1.0, 0.0},
                                         {10, 10, 1}, vel, 1.2, 1.0, 0, 5, 1, initialTemp);
    EXPECT_NEAR(linCelInitTemp.calculateTemperature(), initialTemp, 5.0);
    linCelInitTemp.simulateParticles();
    // after the simulation loop, the system should still have a temperature around 40 degrees
    EXPECT_NEAR(linCelInitTemp.calculateTemperature(), targetTemp, 5.0);
}

//////////////////////////////////////////////////// TESTS FOR THE THERMOSTAT END //////////////////////////////////////



//////////////////////////////////////////////////// TESTS FOR THE RAYLEIGH TAYLOR INSTABILITY BEGIN ///////////////////

//TODO: implement
TEST_F(MolSimTest, testScaleVelocity) {
    
}

// We used the following two tests for print debugging
/*
TEST(MolSimTests, testPeriodicBoundaryToTheLeft) {
    ParticleContainerLinCel linCelTest{0.1, 100, 1, {3.0, 3.0, 1.0}, "pppppp", 1.0};
    linCelTest.add({0.01, 1.5, 0.5}, {-0.1, 0.0, 0.0}, 1, 0);
    std::cout << "Index: " << linCelTest.getCells().at(36).size() << std::endl;
    // linCelTest.getCells().at()
    linCelTest.simulateParticles();
}

TEST(MolSimTests, testPeriodicBoundaryToTheBottom) {
    ParticleContainerLinCel linCelTest{0.1, 100, 1, {3.0, 3.0, 1.0}, "pppppp", 1.0};
    linCelTest.add({0.5, 0.01, 0.5}, {0.0, -0.1, 0.0}, 1, 0);
    std::cout << "Index: " << linCelTest.getCells().at(36).size() << std::endl;
    // linCelTest.getCells().at()
    linCelTest.simulateParticles();
}*/

TEST_F(MolSimTest, testParticleCreation) {
    Particle a{{1.5, 1.5, 0.5}, {0.0, 0.0, 0.0}, 1, 0};
    //check if the default initialization of epsilon and omega
    EXPECT_EQ(a.getSigma(), 1);
    EXPECT_EQ(a.getEpsilon(), 5);
}


//////////////////////////////////////////////////// TESTS FOR THE RAYLEIGH TAYLOR INSTABILITY END /////////////////////