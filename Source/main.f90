!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> In this example, I will show how to convert a CheSS like segmented
!! csr matrix into the NTPoly format, and back again.
PROGRAM ConversionExample
  !! Local Modules
  USE ConversionModule, ONLY : ChessToNTPoly
  USE MPI
  !! NTPoly Modules
  USE ProcessGridModule, ONLY : ConstructProcessGrid, DestructProcessGrid
  USE PSMatrixModule, ONLY : Matrix_ps, CopyMatrix
  !! BigDFT Modules
  USE dictionaries, ONLY : dictionary, dict_free
  USE futile
  USE sparsematrix_highlevel, ONLY: sparse_matrix_metadata_init_from_file, &
       & sparse_matrix_and_matrices_init_from_file_bigdft
  USE sparsematrix_init, ONLY : init_matrix_taskgroups_wrapper
  USE sparsematrix_memory, ONLY : deallocate_sparse_matrix_metadata, &
       & deallocate_matrices, deallocate_sparse_matrix
  USE sparsematrix_types, ONLY : sparse_matrix_metadata, matrices, sparse_matrix
  USE yaml_parse, ONLY : yaml_cl_parse, yaml_cl_parse_null, &
       & yaml_cl_parse_cmd_line, yaml_cl_parse_free, yaml_cl_parse_option
  IMPLICIT NONE
  !! Calculation Parameters
  CHARACTER(LEN=1024) :: metadata_file
  CHARACTER(LEN=1024) :: matrix_file
  CHARACTER(LEN=1024) :: matrix_format
  INTEGER :: num_procs, my_rank
  INTEGER :: ierr, provided
  !! CheSS Matrices
  TYPE(sparse_matrix_metadata) :: metadata
  TYPE(sparse_matrix), DIMENSION(1) :: smat
  TYPE(matrices) :: mat
  !! NTPoly Matrices
  TYPE(Matrix_ps) :: nt_fock, nt_density

  !! Init
  CALL MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, ierr)
  CALL InitParallel
  CALL f_lib_initialize()
  CALL InitParams

  !! Setup the CheSS matrices.
  CALL sparse_matrix_metadata_init_from_file(TRIM(metadata_file), metadata)
  CALL sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, &
       & TRIM(matrix_file), my_rank, num_procs, MPI_COMM_WORLD, smat(1), &
       & mat, init_matmul=.FALSE.)
  CALL init_matrix_taskgroups_wrapper(my_rank, num_procs, MPI_COMM_WORLD, &
       & .TRUE., 1, smat)

  !! Convert CheSS To NTPoly
  CALL ChessToNTPoly(smat(1), mat, nt_fock)

  !! Replace this Copy with whatever solver operation you wish to perform.
  CALL CopyMatrix(nt_fock, nt_density)
  !
  ! !! Now convert back and print.
  ! CALL NTPolyToSeg(NTDensity, Density, start_col, end_col)

  !! Cleanup
  CALL deallocate_sparse_matrix(smat(1))
  CALL deallocate_matrices(mat)
  CALL deallocate_sparse_matrix_metadata(metadata)
  CALL CleanupParallel

CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Read in the input parameters
  SUBROUTINE InitParams()
    !! Local data
    TYPE(yaml_cl_parse) :: parser
    TYPE(dictionary), POINTER :: options

    !! Setup the parser
    parser=yaml_cl_parse_null()
    CALL yaml_cl_parse_option(parser, 'matrix', 'hamiltonian_sparse.txt', &
         & 'Name of the matrix file')
    CALL yaml_cl_parse_option(parser, 'metadata_file', &
         & 'sparsematrix_metadata.dat', &
         & 'The meta data associated with this matrix.')
    CALL yaml_cl_parse_option(parser, 'matrix_format', 'serial_text', &
         & 'Indicates the matrix format.')
    CALL yaml_cl_parse_cmd_line(parser,args=options)
    CALL yaml_cl_parse_free(parser)

    !! Extract
    matrix_file = options//"matrix"
    metadata_file = options//"metadata_file"
    matrix_format = options//"matrix_format"

    !! Punch out
    IF (my_rank .EQ. 0) THEN
       CALL yaml_mapping_open('Input parameters')
       CALL yaml_map('Matrix file',TRIM(matrix_file))
       CALL yaml_map('Metadata file',TRIM(metadata_file))
       CALL yaml_map('Matrix format',TRIM(matrix_format))
       CALL yaml_mapping_close()
    END IF

    !! Cleanup
    CALL dict_free(options)
  END SUBROUTINE InitParams
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Setup MPI/NTPoly
  SUBROUTINE InitParallel
    !! For this driver program
    CALL MPI_Comm_size(MPI_COMM_WORLD, num_procs, ierr)
    CALL MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)

    !! For NTPoly
    CALL ConstructProcessGrid(MPI_COMM_WORLD)
  END SUBROUTINE InitParallel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Tear down MPI/NTPoly
  SUBROUTINE CleanupParallel
    INTEGER :: ierr
    CALL DestructProcessGrid()
    CALL MPI_Finalize(ierr)
  END SUBROUTINE CleanupParallel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END PROGRAM ConversionExample
