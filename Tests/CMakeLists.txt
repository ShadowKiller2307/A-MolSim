include(GoogleTest)

add_executable(
        MolSimTests
        MolSimTests.cpp # file with our tests
)
target_link_libraries(
        MolSimTests
        GTest::gtest_main
)

gtest_discover_tests(MolSimTests)