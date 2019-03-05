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

target_link_libraries(BigPoly ${NTPOLY_LIBRARIES} ${TOOLCHAIN_LIBS}
                      ${BLAS_LIBRARIES} ${MPI_Fortran_LIBRARIES}
                      ${OpenMP_Fortran_LIBRARIES} ${CheSS_LIBRARIES}
                      ${Atlab_LIBRARIES} ${Futile_LIBRARIES} 
                      ${Yaml_LIBRARIES} Database)
target_include_directories(BigPoly PUBLIC ${NTPOLY_INCLUDE_PATH}
    ${Atlab_INCLUDE_PATH} ${Futile_INCLUDE_PATH} 
    ${MPI_INCLUDE_PATH} ${CheSS_INCLUDE_PATH})

install(TARGETS BigPoly DESTINATION bin)
