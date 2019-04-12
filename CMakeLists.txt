################################################################################
## Basic Setup
cmake_minimum_required (VERSION 3.2)
project("BigPoly")
enable_language(Fortran)
enable_testing()

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/Cmake)

################################################################################
## Futile and Chess Dependencies
find_package(Futile REQUIRED)
find_package(CheSS REQUIRED)
find_package(Atlab REQUIRED)
find_package(YAML REQUIRED)

################################################################################
## NTPoly Dependency
find_package(NTPoly REQUIRED)

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
if (NOT ${BLAS_FOUND})
  set(BLAS_LIBRARIES "")
endif()

################################################################################
subdirs(Database)
subdirs(Source)
subdirs(UnitTests)

################################################################################
## Output Locations
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
