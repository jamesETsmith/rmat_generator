cmake_minimum_required(VERSION 3.20)

project(rmat_generator)

set(CMAKE_EXPORT_COMPILE_COMMANDS
    TRUE
    CACHE BOOL "Export compile commands to build directory" FORCE)

#
# Dependencies
#
include(cmake/CPM.cmake)
CPMAddPackage("gh:nlohmann/json@3.10.5")
CPMAddPackage("gh:doctest/doctest@2.4.11")
CPMAddPackage("gh:p-ranav/argparse@3.0")
CPMAddPackage("gh:imneme/pcg-cpp#428802d1a5634f96bcd0705fab379ff0113bcf13")
message(STATUS "HERE ${pcg-cpp_SOURCE_DIR}")

#
# Library
#

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

add_library(rmat_generator INTERFACE)
target_include_directories(rmat_generator INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include ${pcg-cpp_SOURCE_DIR}/include)

#
# Benchmarks and demos
#
add_subdirectory(benchmark)

#
# Test
#
enable_testing()
add_subdirectory(test)

