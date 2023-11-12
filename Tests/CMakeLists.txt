# collect all cpp files
file(GLOB_RECURSE MY_TEST_SRC
    "*.cpp"
    "*.h"
     #[[../src/Particle.h
     ../src/ParticleContainer.h
     ../src/HelperFunctions.h]]
   # ../src/Particle.h
    #../src/ParticleGenerator.h
    #../src/ParticleContainer.h
)

add_executable(MolSimTests ${MY_TEST_SRC}
        ../src/Particle.h
        ../src/ParticleContainer.h
        ../src/HelperFunctions.h
        ../src/Particle.cpp
        ../src/ParticleContainer.cpp
        ../src/HelperFunctions.cpp
        ../src/ForceCalculator.cpp
        ../src/ForceCalculator.h
        ../src/ForceV1.cpp
        ../src/ForceV1.h
        ../src/LennardJonesForce.h
        ../src/LennardJonesForce.cpp
        ../src/FileReader.cpp
        ../src/FileReader.h)

add_test(NAME MolSimTests
         COMMAND MolSimTests)

target_include_directories(MolSimTests PUBLIC ../libs/libxsd PRIVATE ../src)

target_link_libraries(MolSimTests PUBLIC xerces-c GTest::gtest_main PRIVATE spdlog::spdlog nlohmann_json::nlohmann_json)

# activate all compiler warnings. Clean up your code :P
# depending on the compiler different flags are used
target_compile_options(MolSim PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wall> $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>)

gtest_discover_tests(MolSimTests)