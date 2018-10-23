################################################################################
set(sources main.f90 ConversionModule.f90 FakeChessModule.f90)

################################################################################
## Fortran Version
add_executable(convert ${sources})
target_link_libraries(convert ${NTPoly_LIBRARIES} ${Futile_LIBRARIES}
                      ${MPI_Fortran_LIBRARIES} ${OpenMP_Fortran_LIBRARIES}
                      ${BLAS_LIBRARIES} ${Yaml_LIBRARIES} ${CheSS_LIBRARIES})
target_include_directories(convert PUBLIC
  ${NTPoly_INCLUDE_PATH} ${Futile_INCLUDE_PATH} ${MPI_INCLUDE_PATH}
  ${CheSS_INCLUDE_PATH})
install(TARGETS convert DESTINATION bin)
