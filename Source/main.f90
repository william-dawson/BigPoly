!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> In this example, I will show how to convert a CheSS like segmented
!! csr matrix into the NTPoly format, and back again.
PROGRAM ConversionExample
  USE ConversionModule, ONLY : NTPolyToSeg, SegToNTPoly
  !! NTPoly Modules
  USE ProcessGridModule, ONLY : ConstructProcessGrid, DestructProcessGrid
  USE PSMatrixModule, ONLY : Matrix_ps, CopyMatrix, PrintMatrix
  USE TripletListModule, ONLY : TripletList_r
  !! BigDFT Modules
  USE dictionaries, ONLY : dictionary, dict_free
  USE futile
  USE sparsematrix_highlevel, ONLY: sparse_matrix_metadata_init_from_file, &
       & sparse_matrix_and_matrices_init_from_file_bigdft, &
       & ccs_data_from_sparse_matrix
  USE sparsematrix_init, ONLY : init_matrix_taskgroups_wrapper
  USE sparsematrix_memory, ONLY : deallocate_sparse_matrix_metadata, &
       & deallocate_matrices, deallocate_sparse_matrix
  USE sparsematrix_types, ONLY : sparse_matrix_metadata, matrices, sparse_matrix
  USE yaml_parse, ONLY : yaml_cl_parse, yaml_cl_parse_null, &
       & yaml_cl_parse_cmd_line, yaml_cl_parse_free, yaml_cl_parse_option
  USE MPI
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
  !! CCS Matrix
  INTEGER, POINTER, DIMENSION(:) :: row_ind, col_ptr

  !! Init
  CALL MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, ierr)
  CALL f_lib_initialize()
  CALL InitParams
  CALL InitParallel

  !! Setup the CheSS matrices.
  CALL sparse_matrix_metadata_init_from_file(TRIM(metadata_file), metadata)
  CALL sparse_matrix_and_matrices_init_from_file_bigdft(matrix_format, &
       & TRIM(matrix_file), my_rank, num_procs, MPI_COMM_WORLD, smat(1), &
       & mat, init_matmul=.FALSE.)
  CALL init_matrix_taskgroups_wrapper(my_rank, num_procs, MPI_COMM_WORLD, &
       & .TRUE., 1, smat)

  !! Convert to CCS Matrix
  ! col_ptr = f_malloc(smat(1)%nfvctr,id='col_ptr')
  ! row_ind = f_malloc(smat(1)%nvctr,id='row_ind')
  ALLOCATE(col_ptr(smat(1)%nfvctr))
  ALLOCATE(row_ind(smat(1)%nvctr))
  CALL ccs_data_from_sparse_matrix(smat(1), row_ind, col_ptr)

  WRITE(*,*) mat%matrix_compr

  !! The Hamiltonian is a random matrix to start
  ! MatrixSize = options//'matrix_size'
  ! Hamiltonian = SegMat_t(MatrixSize, MatrixSize)
  ! CALL Hamiltonian%FillRandom()
  ! IF (my_rank .EQ. 0) THEN
  !    CALL Hamiltonian%print
  ! END IF
  !
  ! !! The density matrix is empty to start, but here we will construct
  ! !! an empty matrix with a suitable partitioning.
  ! cols_per_proc = MatrixSize/num_procs
  ! start_col = cols_per_proc*(my_rank-1) + 1
  ! end_col = start_col + cols_per_proc - 1
  ! IF (my_rank .EQ. num_procs) THEN
  !    end_col = MatrixSize
  ! END IF
  ! Density = SegMat_t(MatrixSize, cols_per_proc)
  !
  ! !! Now we will demonstrate the conversions. First, convert to NTPoly.
  ! CALL SegToNTPoly(Hamiltonian, NTHamiltonian)
  ! CALL PrintMatrix(NTHamiltonian)
  !
  ! !! Replace this Copy with whatever solver operation you wish to perform.
  ! CALL CopyMatrix(NTHamiltonian, NTDensity)
  !
  ! !! Now convert back and print.
  ! CALL NTPolyToSeg(NTDensity, Density, start_col, end_col)

  !! Cleanup
  ! CALL f_free_ptr(row_ind)
  ! CALL f_free_ptr(col_ptr)
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
    CALL yaml_mapping_open('Input parameters')
    CALL yaml_map('Matrix file',TRIM(matrix_file))
    CALL yaml_map('Metadata file',TRIM(metadata_file))
    CALL yaml_map('Matrix format',TRIM(matrix_format))
    CALL yaml_mapping_close()

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
