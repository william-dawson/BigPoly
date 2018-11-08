################################################################################
set(sources
    main.f90
    AutoYamlArgsModule.f90
    ConversionModule.f90
    MatrixInverseModule.f90
    MultiplyMatricesModule.f90
    ${CMAKE_CURRENT_BINARY_DIR}/Database.f90
)

set (generate
     ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/generatefortran.py
     ${CMAKE_CURRENT_BINARY_DIR}/../Database/input.yaml
     ${CMAKE_CURRENT_BINARY_DIR}/Database.f90
)

add_custom_command(OUTPUT Database.f90
    COMMAND python ${generate}
    DEPENDS ../Database/generatefortran.py ../Database/input.yaml)

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
