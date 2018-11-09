################################################################################
execute_process(COMMAND ${MPIEXEC} --version
                OUTPUT_VARIABLE mpiversion OUTPUT_STRIP_TRAILING_WHITESPACE)
if (${mpiversion} MATCHES "OpenRTE")
  set(oversubscribe "--oversubscribe")
else()
  set(oversubscribe "")
endif()

################################################################################
set(convertargs "--action=convert_matrix_format"
    "--infile=${CMAKE_BINARY_DIR}/UnitTests/hamiltonian_sparse.txt"
    "--outfile=${CMAKE_BINARY_DIR}/UnitTests/output.txt")
add_test(Convert1 mpirun -np 1 ${oversubscribe}
         ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/BigPoly ${convertargs})
add_test(Convert2 mpirun -np 2 ${oversubscribe}
         ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/BigPoly ${convertargs})
add_test(Convert4 mpirun -np 4 ${oversubscribe}
         ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/BigPoly ${convertargs})

set(multiplyargs "--action=multiply_matrices"
    "--infile=${CMAKE_BINARY_DIR}/UnitTests/hamiltonian_sparse.txt"
    "--infile2=${CMAKE_BINARY_DIR}/UnitTests/hamiltonian_sparse.txt"
    "--outfile=${CMAKE_BINARY_DIR}/UnitTests/output.txt")
add_test(Multiply1 mpirun -np 1 ${oversubscribe}
         ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/BigPoly ${multiplyargs})
add_test(Multiply2 mpirun -np 2 ${oversubscribe}
         ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/BigPoly ${multiplyargs})
add_test(Multiply4 mpirun -np 4 ${oversubscribe}
         ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/BigPoly ${multiplyargs})

set(invertargs "--action=matrix_inverse"
    "--infile=${CMAKE_BINARY_DIR}/UnitTests/hamiltonian_sparse.txt"
    "--outfile=${CMAKE_BINARY_DIR}/UnitTests/output.txt")
add_test(Invert1 mpirun -np 1 ${oversubscribe}
         ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/BigPoly ${invertargs})
add_test(Invert2 mpirun -np 2 ${oversubscribe}
         ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/BigPoly ${invertargs})
add_test(Invert3 mpirun -np 4 ${oversubscribe}
         ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/BigPoly ${invertargs})

set(spillageargs "--action=compute_spillage"
   "--infile=${CMAKE_BINARY_DIR}/UnitTests/hamiltonian_sparse.txt"
   "--infile2=${CMAKE_BINARY_DIR}/UnitTests/overlap_sparse.txt"
   "--outfile=${CMAKE_BINARY_DIR}/UnitTests/output.txt"
   "--outfile2=${CMAKE_BINARY_DIR}/UnitTests/output2.txt")
add_test(Spillage1 mpirun -np 1 ${oversubscribe}
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/BigPoly ${spillageargs})
add_test(Spillage2 mpirun -np 2 ${oversubscribe}
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/BigPoly ${spillageargs})
add_test(Spillage4 mpirun -np 4 ${oversubscribe}
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/BigPoly ${spillageargs})

#add_test(Help mpirun -np 1 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/BigPoly --help)

################################################################################
configure_file(hamiltonian_sparse.txt hamiltonian_sparse.txt COPYONLY)
configure_file(overlap_sparse.txt overlap_sparse.txt COPYONLY)
