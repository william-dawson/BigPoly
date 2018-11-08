!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> In this example, I will show how to convert a CheSS like segmented
!! csr matrix into the NTPoly format, and back again.
PROGRAM ConversionExample
  !! Local Modules
  USE ConversionModule, ONLY : ConvertDriver
  USE MatrixPowerModule, ONLY : PowerDriver
  USE MultiplyMatricesModule, ONLY : MultiplyDriver
  !! NTPoly Modules
  USE ProcessGridModule, ONLY : ConstructProcessGrid, DestructProcessGrid
  !! BigDFT Modules
  USE dictionaries
  USE futile
  USE yaml_output, ONLY : yaml_dict_dump_all
  USE yaml_parse, ONLY : yaml_cl_parse, yaml_cl_parse_null, &
       & yaml_parse_from_file, yaml_parse_from_string
  !! ETC
  USE MPI
  IMPLICIT NONE
  !! Calculation Parameters
  INTEGER :: num_procs, my_rank
  INTEGER :: ierr, provided
  TYPE(dictionary), POINTER :: action
  TYPE(dictionary), POINTER :: options

  !! Init
  CALL MPI_Init_thread(MPI_THREAD_SERIALIZED, provided, ierr)
  CALL InitParallel
  CALL f_lib_initialize()
  CALL InitParams

  !! Call A Routine based on the Action Given
  SELECT CASE(TRIM(dict_value(action//'action')))
  CASE("matrixpower")
     ! CALL PowerDriver(options)
  CASE("multiply_matrices")
     CALL MultiplyDriver(options)
  CASE("convert_matrix_format")
     CALL ConvertDriver(options)
  END SELECT

  !! Cleanup
  CALL CleanupParallel
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Read in the input parameters
  SUBROUTINE InitParams()
    !! Local data
    TYPE(yaml_cl_parse) :: parser
    TYPE(dictionary), POINTER :: spec
    TYPE(dictionary), POINTER :: args
    !! The input parameters
    CHARACTER(LEN=*), PARAMETER :: fname = &
         & "/Users/dawson/Desktop/ConversionPoly/Build/Database/input.yaml"
    !! Temporary Variables
    TYPE(dictionary), POINTER :: dtmp, it

    !! Setup the initial spec from file.
    CALL yaml_parse_from_file(dtmp, fname)
    spec => dtmp .pop. 0
    CALL dict_free(dtmp)

    !! The first input option is what action to take.
    parser=yaml_cl_parse_null()
    it => dict_iter(spec)
    CALL yaml_cl_parse_option(parser, "action", TRIM(dict_key(it)), &
         'What action to take')
    CALL yaml_cl_parse_cmd_line(parser,args=action)
    CALL yaml_cl_parse_free(parser)
    CALL yaml_map('Parsed action',action)

    !! Next we parse the options associated with that action.
    parser=yaml_cl_parse_null()
    args => (spec//dict_value(action//'action'))//'args'
    it => dict_iter(args)
    DO WHILE(ASSOCIATED(it))
       CALL yaml_cl_parse_option(parser, TRIM(dict_key(it)), &
            & TRIM(dict_value(it//'default')), &
            & TRIM(dict_value(it//'shorthelp')))
       it => dict_next(it)
    END DO
    CALL yaml_cl_parse_cmd_line(parser,args=options)
    CALL yaml_cl_parse_free(parser)
    CALL yaml_map('Parsed options',options)
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
