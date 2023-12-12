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

# See https://best.openssf.org/Compiler-Hardening-Guides/Compiler-Options-Hardening-Guide-for-C-and-C++
# -fstrict-flex-arrays=3 requires gcc 13+ and llvm 16+
# -D_GLIBCXX_ASSERTIONS -D_LIBCPP_ASSERT 
# -mbranch-protection=standard
set(COMPILE_TIME_HARDENING_FLAGS "-Wall -Wformat=2 -Wconversion -Wtrampolines -Wimplicit-fallthrough")
set(RUN_TIME_HARDENING_FLAGS "-U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=3 -D_GLIBCXX_ASSERTIONS -fstack-clash-protection -fcf-protection=full ")

set(CMAKE_CXX_FLAGS "-g -O3 -march=native ${COMPILE_TIME_HARDENING_FLAGS} ${RUN_TIME_HARDENING_FLAGS}")

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

