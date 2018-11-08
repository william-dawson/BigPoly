!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This modules has the routines to convert an arbitrary Inverse of a CheSS
!! matrix.
MODULE MatrixInverseModule
  !! Local modules
  USE ConversionModule, ONLY : ChessToNTPoly
  !! BigDFT Modules
  USE dictionaries
  USE futile
  !! NTPoly Modules
  USE DataTypesModule, ONLY : NTREAL
  USE InverseSolversModule, ONLY : Invert
  USE PSMatrixModule, ONLY : Matrix_ps, DestructMatrix, &
       & WriteMatrixToMatrixMarket
  USE SolverParametersModule, ONLY : SolverParameters_t, &
       & DestructSolverParameters
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: InverseDriver
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE InverseDriver(options)
    !> The dictionary specifying the parameters of the calculation.
    TYPE(dictionary), INTENT(IN), POINTER :: options
    !! Local variables
    TYPE(Matrix_ps) :: mat, invmat
    REAL(f_double) :: threshold, convergence
    TYPE(SolverParameters_t) :: param

    !! Read in the input
    CALL ChessToNTPoly(dict_value(options//'matrix_format'), &
         & dict_value(options//'infile'), mat)
    threshold = options//'threshold'
    convergence = options//'convergence_threshold'

    !! Do the inversion
    param = SolverParameters_t()
    param%threshold = REAL(threshold,NTREAL)
    param%converge_diff = REAL(convergence,NTREAL)
    param%be_verbose = .TRUE.
    CALL Invert(mat, invmat, param)

    !! Write To File
    CALL WriteMatrixToMatrixMarket(invmat, dict_value(options//'outfile'))

    !! Cleanup
    CALL DestructSolverParameters(param)
    CALL DestructMatrix(mat)
    CALL DestructMatrix(invmat)
  END SUBROUTINE InverseDriver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE MatrixInverseModule
