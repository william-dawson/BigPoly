################################################################################
## Basic Setup
cmake_minimum_required (VERSION 3.2)
project("Conversion Poly")
enable_language(Fortran)
enable_testing()

################################################################################
## NTPoly, Futile, and Chess Dependencies
find_library(NTPoly_LIBRARIES NTPoly)
find_path(NTPoly_INCLUDE_PATH "psmatrixmodule.mod")
find_library(Futile_LIBRARIES futile-1)
find_path(Futile_INCLUDE_PATH "dictionaries.mod")

################################################################################
## Packages
find_package(MPI REQUIRED)
find_package(OpenMP)
find_package(BLAS)

################################################################################
subdirs(Source)

################################################################################
## Output Locations
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

################################################################################
add_test(Test1 mpirun -np 1 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/example_fortran)
