cmake_minimum_required(VERSION 3.16)

# define project name, version
project(PSEMolDyn_GroupA)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# set Release as the default build type if it is not yet set.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

file(MAKE_DIRECTORY "output")

add_subdirectory(src)
add_subdirectory(Tests EXCLUDE_FROM_ALL)

find_package(Doxygen)

#include(doxygen)
# TODO: Include doxygen 