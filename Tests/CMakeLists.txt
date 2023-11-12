# collect all cpp files
file(GLOB_RECURSE MY_SRC
    "*.cpp"
    "*.h"
    "../src/*.cpp"
    "../src/*.h"
)

add_executable(MolSimTests ${MY_SRC})

target_include_directories(MolSimTests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs/libxsd PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# according to http://google.github.io/googletest/quickstart-cmake.html
include(FetchContent)
FetchContent_Declare(googletest URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)
enable_testing()

target_link_libraries(MolSimTests PUBLIC xerces-c GTest::gtest_main PRIVATE nlohmann_json::nlohmann_json)

# activate all compiler warnings. Clean up your code :P
# depending on the compiler different flags are used
target_compile_options(MolSim PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wall> $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>)

gtest_discover_tests(MolSimTests)