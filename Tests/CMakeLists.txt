# collect all cpp files
file(GLOB_RECURSE MY_TEST_SRC
    "*.cpp"
    "*.h"
    "../src/*.cpp"
    "../src/*.h"
)

add_executable(MolSimTests ${MY_TEST_SRC})

target_include_directories(MolSimTests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_link_libraries(MolSimTests PUBLIC xerces-c GTest::gtest_main PRIVATE nlohmann_json::nlohmann_json)

# activate all compiler warnings. Clean up your code :P
# depending on the compiler different flags are used
target_compile_options(MolSim PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wall> $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>)

gtest_discover_tests(MolSimTests)