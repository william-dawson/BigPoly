################################################################################
## Basic Setup
cmake_minimum_required (VERSION 3.2)
project("Conversion Poly")
enable_language(Fortran)
enable_testing()

################################################################################
## Futile and Chess Dependencies
find_library(Futile_LIBRARIES futile-1)
find_path(Futile_INCLUDE_PATH "dictionaries.mod")
find_library(Yaml_LIBRARIES yaml)
find_library(CheSS_LIBRARIES CheSS-1)
find_path(CheSS_INCLUDE_PATH "sparsematrix.mod")

################################################################################
## NTPoly Dependency
find_library(NTPOLY_LIBRARIES NTPoly)
find_path(NTPOLY_INCLUDE_PATH "psmatrixmodule.mod")
MESSAGE(${NTPOLY_LIBRARIES})

################################################################################
## Output Locations
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)

################################################################################
## Packages
find_package(MPI REQUIRED)
find_package(OpenMP)
find_package(BLAS)

################################################################################
subdirs(Database)
subdirs(Source)
subdirs(UnitTests)

################################################################################
## Output Locations
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
