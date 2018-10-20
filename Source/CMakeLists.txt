################################################################################
set(sources main.f90 ConversionModule.f90 FakeChessModule.f90)

################################################################################
## Fortran Version
add_executable(example_fortran ${sources})
target_link_libraries(example_fortran ${NTPoly_LIBRARIES} ${Futile_LIBRARIES}
                      ${MPI_Fortran_LIBRARIES} ${OpenMP_Fortran_LIBRARIES}
                      ${BLAS_LIBRARIES})
target_include_directories(example_fortran PUBLIC
  ${NTPoly_INCLUDE_PATH} ${Futile_INCLUDE_PATH} ${MPI_INCLUDE_PATH})
install(TARGETS example_fortran DESTINATION bin)
