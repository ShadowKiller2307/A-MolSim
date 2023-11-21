# collect all cpp files
file(GLOB_RECURSE MY_SRC
    "*.cpp"
    "*.h"
)

# create make target
add_executable(MolSim ${MY_SRC}
        ParticleContainer.cpp
        ParticleContainer.h
        ParticleContainerLC.cpp
        ParticleContainerLC.h
        BoundaryCondition.cpp
        BoundaryCondition.h
)

target_include_directories(MolSim PUBLIC ../libs/libxsd PRIVATE .)

target_link_libraries(MolSim PUBLIC xerces-c PRIVATE spdlog::spdlog nlohmann_json::nlohmann_json)

# activate all compiler warnings. Clean up your code :P
# depending on the compiler different flags are used
target_compile_options(MolSim PRIVATE $<$<CXX_COMPILER_ID:GNU>:-Wall> $<$<CXX_COMPILER_ID:Intel>:-w3 -wd383,981,1418,1572,2259>)