!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> In this example, I will show how to convert a CheSS like segmented
!! csr matrix into the NTPoly format, and back again.
PROGRAM ConversionExample
  !! Local Modules
  USE AutoYamlArgsModule, ONLY : GetArgs
  USE ConversionModule, ONLY : ConvertDriver
  USE MatrixInverseModule, ONLY : InverseDriver
  USE MultiplyMatricesModule, ONLY : MultiplyDriver
  !! NTPoly Modules
  USE ProcessGridModule, ONLY : ConstructProcessGrid, DestructProcessGrid
  !! BigDFT Modules
  USE dictionaries
  USE futile
  !! ETC
  USE MPI
  IMPLICIT NONE
  !! Calculation Parameters
  INTEGER :: num_procs, my_rank
  INTEGER :: ierr, provided
  TYPE(dictionary), POINTER :: options

  !! Init
  CALL MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, ierr)
  CALL InitParallel
  CALL f_lib_initialize()
  CALL GetArgs(options, my_rank)

  !! Call A Routine based on the Action Given
  SELECT CASE(TRIM(dict_value(options//'action')))
  CASE("matrix_inverse")
     CALL InverseDriver(options)
  CASE("multiply_matrices")
     CALL MultiplyDriver(options)
  CASE("convert_matrix_format")
     CALL ConvertDriver(options)
  END SELECT

  !! Cleanup
  CALL dict_free(options)
  CALL CleanupParallel
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
