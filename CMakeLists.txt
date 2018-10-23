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
find_library(Yaml_LIBRARIES yaml)
find_library(CheSS_LIBRARIES CheSS-1)
find_path(CheSS_INCLUDE_PATH "sparsematrix.mod")

################################################################################
## Packages
find_package(MPI REQUIRED)
find_package(OpenMP)
find_package(BLAS)

################################################################################
subdirs(Source)
subdirs(UnitTests)

################################################################################
## Output Locations
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
