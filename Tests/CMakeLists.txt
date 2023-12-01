# collect all cpp files
file(GLOB_RECURSE MY_TEST_SRC RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    "*.cpp"
    "*.h"
    "../src/*.cpp"
    "../src/*.h"
)

#remove the MolSim.cpp and MolSim.h as they are not needed for the tests
string(REPLACE "../src/MolSim.cpp" "" MY_TEST_SRC "${MY_TEST_SRC}")
string(REPLACE "../src/MolSim.h" "" MY_TEST_SRC "${MY_TEST_SRC}")

add_executable(MolSimTests ${MY_TEST_SRC}
        ../src/xmlSchema/schema.cpp)

target_include_directories(MolSimTests PUBLIC ../libs/libxsd PRIVATE ../src)

target_link_libraries(MolSimTests PUBLIC xerces-c GTest::gtest_main PRIVATE spdlog::spdlog nlohmann_json::nlohmann_json)

# activate all compiler warnings. Clean up your code :P
# depending on the compiler different flags are used
target_compile_options(MolSim PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wall> $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>)

gtest_discover_tests(MolSimTests)
gtest_add_tests(TARGET MolSimTests
                TEST_SUFFIX .noArgs
                TEST_LIST   noArgsTests
)
add_test(NAME MolSimTests
        COMMAND $<TARGET_FILE:MolSimTests>)