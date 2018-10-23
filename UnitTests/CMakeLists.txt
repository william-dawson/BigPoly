
################################################################################
execute_process(COMMAND ${MPIEXEC} --version
                OUTPUT_VARIABLE mpiversion OUTPUT_STRIP_TRAILING_WHITESPACE)
if (${mpiversion} MATCHES "OpenRTE")
  set(oversubscribe "--oversubscribe")
else()
  set(oversubscribe "")
endif()

################################################################################
set(testargs "--matrix
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/UnitTests/hamiltonian_sparse.txt
  --metadata
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/UnitTests/sparsematrix_metadata.dat")

add_test(Test1 mpirun -np 1 ${oversubscribe}
         ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/convert)
add_test(Test2 mpirun -np 2 ${oversubscribe}
         ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/convert)
add_test(Test4 mpirun -np 4 ${oversubscribe}
         ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/convert)

################################################################################
configure_file(hamiltonian_sparse.txt hamiltonian_sparse.txt COPYONLY)
configure_file(sparsematrix_metadata.dat sparsematrix_metadata.dat COPYONLY)