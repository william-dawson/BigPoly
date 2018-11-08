################################################################################
# Required Packages
find_package(PythonInterp REQUIRED)

################################################################################
# Files for doing the code generation
configure_file(input.yaml input.yaml COPYONLY)
configure_file(generatefortran.py
               ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/generatefortran.py COPYONLY)

################################################################################
# Actual code generation
set (generate
    ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/generatefortran.py
    ${CMAKE_CURRENT_BINARY_DIR}/input.yaml
    ${CMAKE_CURRENT_BINARY_DIR}/Database.f90
)

add_custom_command(OUTPUT Database.f90
   COMMAND ${PYTHON_EXECUTABLE} ${generate}
   DEPENDS generatefortran.py input.yaml)

################################################################################
# Build the library
add_library(Database STATIC Database.f90)
include(GNUInstallDirs)
install(TARGETS Database
        EXPORT Database-export
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
