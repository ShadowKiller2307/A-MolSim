cmake_minimum_required(VERSION 3.16)

# define project name, version
project(PSEMolDyn_GroupA)

# let ccmake and cmake-gui offer the default build type options
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug;Release;RelWithDebInfo;MinSizeRel")

# set Release as the default build type if it is not yet set.
#if(NOT CMAKE_BUILD_TYPE)
#    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
#            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
#endif(NOT CMAKE_BUILD_TYPE)
#set(CMAKE_BUILD_TYPE RelWithDebInfo)
set(CMAKE_BUILD_TYPE Debug)

file(MAKE_DIRECTORY "output")
file(MAKE_DIRECTORY "logs")

include(FetchContent)
include(CTest)

# according to http://google.github.io/googletest/quickstart-cmake.html
FetchContent_Declare(googletest URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip)
FetchContent_MakeAvailable(googletest)
include(CTest)
include(GoogleTest)
enable_testing()
FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz)
FetchContent_MakeAvailable(json)
FetchContent_Declare(spdlog GIT_REPOSITORY https://github.com/gabime/spdlog.git GIT_TAG v1.x)
FetchContent_MakeAvailable(spdlog)

add_subdirectory(src)
add_subdirectory(Tests EXCLUDE_FROM_ALL)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

option(BUILD_DOXY "Create Doxygen Documentation" ON)

if (BUILD_DOXY)
    find_package(Doxygen)

    if (DOXYGEN_FOUND)
        add_custom_target(
                doc_doxygen
                COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
        set_target_properties(doc_doxygen PROPERTIES EXCLUDE_FROM_ALL TRUE)
    endif (DOXYGEN_FOUND)
endif (BUILD_DOXY)