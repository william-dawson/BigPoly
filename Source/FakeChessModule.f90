!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This module contains fake data structures imitating CheSS
MODULE FakeChessModule
  IMPLICIT NONE
  !> An example implementation of the segmented csr matrix.
  TYPE, PUBLIC :: SegMat_t
     !> The number of columns in a matrix
     INTEGER :: num_columns
     !> The number of rows in a matrix
     INTEGER :: num_rows
     !> How many segments are in each column.
     INTEGER, DIMENSION(:), ALLOCATABLE :: outer_index
     !> What the starting index is if each segment in a column.
     INTEGER, DIMENSION(:), ALLOCATABLE :: segment_index
     !> The length of each segement.
     INTEGER, DIMENSION(:), ALLOCATABLE :: segment_len
     !> The actual data.
     REAL(8), DIMENSION(:), ALLOCATABLE :: values
   CONTAINS
     PROCEDURE :: FillRandom
     PROCEDURE :: PRINT => PrintMat
     FINAL :: Final_segmat
  END TYPE SegMat_t
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTERFACE SegMat_t
     MODULE PROCEDURE segmat_init
  END INTERFACE
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Initializes a CSR matrix with a random matrix. You should ignore this
  !! code as it's just an example.
  FUNCTION segmat_init(rows, columns) RESULT(this)
    TYPE(SegMat_t) :: this
    INTEGER, INTENT(IN) :: rows
    INTEGER, INTENT(IN) :: columns

    !! Process Args
    this%num_rows = rows
    this%num_columns = columns
  END FUNCTION segmat_init
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out a matrix to the console.
  !! (As a sanity check)
  SUBROUTINE PrintMat(this)
    !> The matrix to print out.
    CLASS(SegMat_t), INTENT(IN) :: this
    !! Local Variables
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: densemat
    INTEGER :: OO, SS
    INTEGER :: segstart, segend

    !! We'll build a dense matrix to print
    CALL CsrToDense(this, densemat)

    CALL PrintDense(densemat)

    !! Cleanup
    DEALLOCATE(densemat)

  END SUBROUTINE PrintMat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Fill with random entries.
  SUBROUTINE FillRandom(this)
    !> The matrix to fill.
    CLASS(SegMat_t), INTENT(INOUT) :: this
    !! Local Variables
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: densemat
    INTEGER :: II
    REAL(8) :: rand_temp

    CALL Final_segmat(this)

    !! Dense random matrix.
    CALL ComputeRandomMask(densemat, this%num_rows, this%num_columns)
    CALL DenseToCsr(densemat, this)

    !! Next we compute the actual values.
    DO II = 1, SIZE(this%values)
       CALL RANDOM_NUMBER(rand_temp)
       this%values(II) = rand_temp
    END DO

    !! Cleanup
    DEALLOCATE(densemat)
  END SUBROUTINE FillRandom
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Destructor for the CSR Matrix
  SUBROUTINE Final_segmat(this)
    TYPE(SegMat_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%outer_index)) DEALLOCATE(this%outer_index)
    IF (ALLOCATED(this%segment_index)) DEALLOCATE(this%segment_index)
    IF (ALLOCATED(this%segment_len)) DEALLOCATE(this%segment_len)
    IF (ALLOCATED(this%values)) DEALLOCATE(this%values)
  END SUBROUTINE Final_segmat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE CsrToDense(segmat, dmat)
    !> The matrix to print out.
    TYPE(SegMat_t), INTENT(IN) :: segmat
    REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: dmat
    !! Local variables
    INTEGER :: OO, SS, VV
    INTEGER :: segstart, segend
    INTEGER :: iptr

    !! We'll build a dense matrix to print
    ALLOCATE(dmat(segmat%num_rows, segmat%num_columns))
    dmat = 0

    !! Loop over outer index
    iptr = 1
    DO OO = 1, segmat%num_columns
       !! Loop over segments
       DO SS = segmat%outer_index(OO), segmat%outer_index(OO+1)-1
          segstart = segmat%segment_index(SS)
          DO VV = 1, segmat%segment_len(SS)
             dmat(segstart+VV-1,OO) = segmat%values(iptr)
             iptr = iptr + 1
          END DO
       END DO
    END DO
  END SUBROUTINE CsrToDense
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert a dense matrix to a csr matrix
  SUBROUTINE DenseToCsr(dmat, segmat)
    !> Dense matrix to convert
    REAL(8), DIMENSION(:,:), INTENT(IN) :: dmat
    !> CSR matri to convert to
    TYPE(SegMat_t) :: segmat
    !! Local variables
    INTEGER :: II, JJ
    INTEGER, DIMENSION(:), ALLOCATABLE :: segment_end
    INTEGER :: nnz
    INTEGER :: eptr, sptr

    segmat = SegMat_t(SIZE(dmat,DIM=1), SIZE(dmat,DIM=2))

    !! Compute outer index.
    ALLOCATE(segmat%outer_index(segmat%num_columns+1))
    segmat%outer_index = 0
    segmat%outer_index(1) = 1
    DO JJ = 2, segmat%num_columns+1
       segmat%outer_index(JJ) = segmat%outer_index(JJ-1)
       DO II = 1, segmat%num_rows
          IF (dmat(II,JJ-1) .EQ. 1) THEN
             IF (II .EQ. 1) THEN
                segmat%outer_index(JJ) = segmat%outer_index(JJ) + 1
             ELSE IF (dmat(II-1,JJ-1) .EQ. 0) THEN
                segmat%outer_index(JJ) = segmat%outer_index(JJ) + 1
             END IF
          END IF
       END DO
    END DO
    nnz = segmat%outer_index(segmat%num_columns+1) - 1

    !! Count the segment lengths and indices
    ALLOCATE(segmat%segment_len(nnz))
    segmat%segment_len = 0
    ALLOCATE(segmat%segment_index(nnz))
    segmat%segment_index = 0
    ALLOCATE(segment_end(nnz))
    sptr = 1
    eptr = 1
    DO JJ = 1, segmat%num_columns
       DO II = 1, segmat%num_rows
          IF (dmat(II,JJ) .EQ. 1) THEN
             IF (segmat%num_rows .EQ. 1) THEN
                segmat%segment_index(sptr) = 1
                segment_end(eptr) = 1
                sptr = sptr + 1
                eptr = eptr + 1
             ELSE
                !! Check Segment Start
                IF (II .EQ. 1) THEN
                   segmat%segment_index(sptr) = II
                   sptr = sptr + 1
                ELSE IF (dmat(II-1,JJ) .EQ. 0) THEN
                   segmat%segment_index(sptr) = II
                   sptr = sptr + 1
                END IF
                !! Check Segment End
                IF (II .EQ. segmat%num_rows) THEN
                   segment_end(eptr) = II
                   eptr = eptr + 1
                ELSE IF (dmat(II+1,JJ) .EQ. 0) THEN
                   segment_end(eptr) = II
                   eptr = eptr + 1
                END IF
             END IF
          END IF
       END DO
    END DO

    DO II = 1, nnz
       segmat%segment_len(II) = segment_end(II) - segmat%segment_index(II) + 1
    END DO

    !! Allocate Values
    ALLOCATE(segmat%values(SUM(segmat%segment_len)))
    segmat%values = 1

    !! Cleanup
    DEALLOCATE(segment_end)

  END SUBROUTINE DenseToCsr
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Computes a random matrix pattern.
  SUBROUTINE ComputeRandomMask(densemat, rows, cols)
    !> Matrix to compute.
    REAL(8), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: densemat
    !> Number of rows of the matrix.
    INTEGER :: rows
    !> Number of columns of the matrix.
    INTEGER :: cols
    REAL(8) :: rand_t
    INTEGER :: II, JJ

    IF (ALLOCATED(densemat)) THEN
       DEALLOCATE(densemat)
    END IF
    ALLOCATE(densemat(rows, cols))
    densemat = 0

    DO II = 1, rows
       DO JJ = 1, cols
          CALL RANDOM_NUMBER(rand_t)
          IF (rand_t > 0.5) THEN
             densemat(II,JJ) = 1
          END IF
       END DO
    END DO
  END SUBROUTINE ComputeRandomMask
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Print out a dense matrix.
  SUBROUTINE PrintDense(mat)
    !> Matrix to print
    REAL(8), DIMENSION(:,:), INTENT(IN) :: mat
    INTEGER :: II, JJ

    DO II = 1, SIZE(mat,DIM=1)
       WRITE(*,*) mat(II,:)
    END DO
  END SUBROUTINE PrintDense
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE FakeChessModule
