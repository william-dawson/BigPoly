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

if (NTPOLY_LIB)
  target_link_libraries(BigPoly NTPoly::NTPoly ${Futile_LIBRARIES}
                        ${MPI_Fortran_LIBRARIES} ${OpenMP_Fortran_LIBRARIES}
                        ${BLAS_LIBRARIES} ${Yaml_LIBRARIES} ${CheSS_LIBRARIES}
                        Database)
  target_include_directories(BigPoly PUBLIC NTPoly::NTPoly
      ${Futile_INCLUDE_PATH} ${MPI_INCLUDE_PATH} ${CheSS_INCLUDE_PATH})
else()
  target_link_libraries(BigPoly ${Futile_LIBRARIES} NTPoly
                        ${MPI_Fortran_LIBRARIES} ${OpenMP_Fortran_LIBRARIES}
                        ${BLAS_LIBRARIES} ${Yaml_LIBRARIES} ${CheSS_LIBRARIES}
                        Database)
  target_include_directories(BigPoly PUBLIC
    ${Futile_INCLUDE_PATH} ${MPI_INCLUDE_PATH}
    ${CheSS_INCLUDE_PATH})
endif()

install(TARGETS BigPoly DESTINATION bin)
