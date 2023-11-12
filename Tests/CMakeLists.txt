# according to http://google.github.io/googletest/quickstart-cmake.html
include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)
enable_testing()



# collect all cpp files
file(GLOB_RECURSE MY_SRC
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
        # header don't need to be included but this might be necessary for some IDEs
        "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
)

set(MOLSIM_SOURCES
        ../src/ParticleContainer.cpp
        ../src/Cuboid.cpp
        ../src/LennardJonesForce.cpp
        ../src/Particle.cpp
        ../src/HelperFunctions.cpp
        ../src/ForceV1.cpp
        ../src/ForceCalculator.cpp
        ../src/FileReader.cpp
        MolSimTests.cpp
)
set(MOLSIM_HEADERS

        ../src/ParticleContainer.h
        ../src/Particle.h
        ../src/LennardJonesForce.h
        ../src/ForceV1.h
        ../src/HelperFunctions.h
        ../src/ForceCalculator.h
        ../src/FileReader.h

)

add_executable(
        MolSimTests ${MY_SRC}
        ${MOLSIM_SOURCES}
        ${MOLSIM_HEADERS}

         # file with our tests
)
target_link_libraries(
        MolSimTests
        GTest::gtest_main
)

target_include_directories(MolSimTests
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

gtest_discover_tests(MolSimTests)