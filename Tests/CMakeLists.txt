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
        "*.cpp"
        # header don't need to be included but this might be necessary for some IDEs
        "*.h"
        "../src/*.cpp"
        "../src/*.h"
)

add_executable(
        MolSimTests ${MY_SRC}
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