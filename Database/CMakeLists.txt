################################################################################
configure_file(input.yaml input.yaml COPYONLY)
configure_file(generatefortran.py
               ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/generatefortran.py COPYONLY)


               set (generate
                    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/generatefortran.py
                    ${CMAKE_CURRENT_BINARY_DIR}/../Database/input.yaml
                    ${CMAKE_CURRENT_BINARY_DIR}/Database.f90
               )

               add_custom_command(OUTPUT Database.f90
                   COMMAND python ${generate}
                   DEPENDS ../Database/generatefortran.py ../Database/input.yaml)
