!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This modules has the routines to multiply two CheSS matrices.
MODULE MultiplyMatricesModule
  !! BigDFT Modules
  USE dictionaries
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: MultiplyDriver
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE MultiplyDriver(options)
    !> The dictionary specifying the parameters of the calculation.
    TYPE(dictionary), INTENT(IN), POINTER :: options
  END SUBROUTINE MultiplyDriver
END MODULE MultiplyMatricesModule
