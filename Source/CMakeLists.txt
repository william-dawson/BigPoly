################################################################################
set(sources
    main.f90
    AutoYamlArgsModule.f90
    ConversionModule.f90
    MatrixInverseModule.f90
    MultiplyMatricesModule.f90
    SpillageModule.f90
)

################################################################################
## Build the executable
add_executable(BigPoly ${sources})

target_link_libraries(BigPoly ${NTPoly_LIBRARIES} ${TOOLCHAIN_LIBS}
                      ${BLAS_LIBRARIES} ${MPI_Fortran_LIBRARIES}
                      ${OpenMP_Fortran_LIBRARIES} ${CheSS_LIBRARIES}
                      ${Atlab_LIBRARIES} ${Futile_LIBRARIES}
                      ${YAML_LIBRARIES} Database)

target_include_directories(BigPoly PUBLIC ${NTPoly_INCLUDE_DIRS}
    ${Atlab_INCLUDE_DIRS} ${Futile_INCLUDE_DIRS}
  ${MPI_INCLUDE_DIRS} ${CheSS_INCLUDE_DIRS} ${MPI_Fortran_INCLUDE_DIRS})

install(TARGETS BigPoly DESTINATION bin)
