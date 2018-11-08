!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This modules has the routines to multiply two CheSS matrices.
MODULE MultiplyMatricesModule
  !! Local modules
  USE ConversionModule, ONLY : ChessToNTPoly
  !! BigDFT Modules
  USE dictionaries
  USE futile
  !! NTPoly Modules
  USE DataTypesModule, ONLY : NTREAL
  USE PSMatrixModule, ONLY : Matrix_ps, DestructMatrix, &
       & WriteMatrixToMatrixMarket
  USE PSMatrixAlgebraModule, ONLY : MatrixMultiply
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: MultiplyDriver
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE MultiplyDriver(options)
    !> The dictionary specifying the parameters of the calculation.
    TYPE(dictionary), INTENT(IN), POINTER :: options
    !! Local variables
    TYPE(Matrix_ps) :: mata, matb, matc
    REAL(f_double) :: threshold

    !! Read in the input
    CALL ChessToNTPoly(dict_value(options//'matrix_format'), &
         & dict_value(options//'infile'), mata)
    CALL ChessToNTPoly(dict_value(options//'matrix_format'), &
         & dict_value(options//'infile2'), matb)
    threshold = options//'threshold'

    !! Do the actual multiplication
    CALL MatrixMultiply(mata, matb, matc, threshold_in=REAL(threshold,NTREAL))

    !! Write To File
    CALL WriteMatrixToMatrixMarket(matc, dict_value(options//'outfile'))

    !! Cleanup
    CALL DestructMatrix(mata)
    CALL DestructMatrix(matb)
    CALL DestructMatrix(matc)
  END SUBROUTINE MultiplyDriver
END MODULE MultiplyMatricesModule
