!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This modules has the routines to compute the matrices needed for the
!! spillage analysis.
MODULE SpillageModule
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
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply
  USE SolverParametersModule, ONLY : SolverParameters_t, &
       & DestructSolverParameters
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: SpillageDriver
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE SpillageDriver(options)
    !> The dictionary specifying the parameters of the calculation.
    TYPE(dictionary), INTENT(IN), POINTER :: options
    !! Local variables
    TYPE(Matrix_ps) :: hmat, smat, sinv, sinvh, sinvh2
    REAL(f_double) :: threshold, convergence
    TYPE(SolverParameters_t) :: param

    !! Read in the input
    CALL ChessToNTPoly(dict_value(options//'matrix_format'), &
         & dict_value(options//'infile'), hmat)
    CALL ChessToNTPoly(dict_value(options//'matrix_format'), &
         & dict_value(options//'infile2'), smat)
    threshold = options//'threshold'
    convergence = options//'convergence_threshold'

    !! Do the inversion
    param = SolverParameters_t()
    param%threshold = REAL(threshold,NTREAL)
    param%converge_diff = REAL(convergence,NTREAL)
    param%be_verbose = .TRUE.
    CALL Invert(smat, sinv, param)

    !! Do the multiplication
    CALL MatrixMultiply(sinv, hmat, sinvh, threshold_in=REAL(threshold,NTREAL))
    CALL MatrixMultiply(sinvh, sinvh, sinvh2, &
         & threshold_in=REAL(threshold,NTREAL))

    !! Write To File
    CALL WriteMatrixToMatrixMarket(sinvh, dict_value(options//'outfile'))
    CALL WriteMatrixToMatrixMarket(sinvh2, dict_value(options//'outfile2'))

    !! Cleanup
    CALL DestructMatrix(hmat)
    CALL DestructMatrix(smat)
    CALL DestructMatrix(sinv)
    CALL DestructMatrix(sinvh)
    CALL DestructMatrix(sinvh2)
  END SUBROUTINE SpillageDriver
END MODULE SpillageModule
