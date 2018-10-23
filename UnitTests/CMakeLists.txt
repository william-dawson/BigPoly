################################################################################
set(testargs "--matrix
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/UnitTests/hamiltonian_sparse.txt
  --metadata
  ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/UnitTests/sparsematrix_metadata.dat")
add_test(Test1 mpirun -np 1 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/convert)
add_test(Test2 mpirun -np 2 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/convert)

################################################################################
configure_file(hamiltonian_sparse.txt hamiltonian_sparse.txt COPYONLY)
configure_file(sparsematrix_metadata.dat sparsematrix_metadata.dat COPYONLY)
