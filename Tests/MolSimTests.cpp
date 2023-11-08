/**
 * @brief Class for Unit Tests which will check the correctness of our code
 */

#include <gtest/gtest.h>
#include "../src/ParticleContainer.h"
#include "../src/HelperFunctions.h"
#include "../src/Particle.h"

// Difference ASSERT vs EXPECT macros
//ASSERT -> fatal failures
// EXPECT -> nonfatal failures

//
// unit Test for the particle container
TEST(ParticleContainerTest, testGetParticles) {
    //ParticleContainer container;  //Hier also noch  undefined reference to particle container aus irgendeinem Grund
    EXPECT_EQ(0, 0); //Normale Tests funktionieren
    //TODO: Error undefined Reference to HelperFunction
    // Es kann noch in der Testklasse nicht auf unsere Klassen aus dem src Ordner zugegriffen werden
//    auto testArray{std::array<double, 3>{0.0, 0.0, 0.0}};
//    HelperFunctions::scalarOperations(testArray, 1, false);
//    EXPECT_EQ(0.0, testArray[0]);
}
