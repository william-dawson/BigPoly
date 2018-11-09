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
find_package(NTPoly QUIET)
IF (NOT NTPOLY_LIB)
  include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
  ExternalProject_Add(NTPoly
    GIT_REPOSITORY    https://github.com/william-dawson/NTPoly.git
    GIT_TAG           origin/CMake
    CMAKE_ARGS        -DFORTRAN_ONLY=Yes
  )
  ExternalProject_Get_Property(NTPoly source_dir)
  find_package(NTPoly HINTS ${source_dir})
ENDIF()

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
