!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This modules has the routines to convert a CSR matrix to NTPoly and
!! back again.
MODULE ConversionModule
  USE FakeChessModule, ONLY : CSRMat_t
  USE PSMatrixModule, ONLY : Matrix_ps, GetMatrixTripletList, &
       & ConstructEmptyMatrix, FillMatrixFromTripletList, GetMatrixBlock
  USE TripletListModule, ONLY : TripletList_r, ConstructTripletList, &
       & AppendToTripletList, DestructTripletList
  USE TripletModule, ONLY : Triplet_r
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: CSRToNTPoly
  PUBLIC :: NTPolyToCSR
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert a CSRMat_t type to NTPoly type
  SUBROUTINE CSRToNTPoly(csrmat, ntpolymat)
    !> The matrix to convert.
    TYPE(CSRMat_t), INTENT(IN) :: csrmat
    !> The output matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: ntpolymat
    !! Local Data
    TYPE(CSRMat_t) :: csrsplit
    TYPE(TripletList_r) :: triplet_list
    TYPE(Triplet_r) :: trip
    INTEGER :: II, JJ, KK
    INTEGER :: iptr
    INTEGER :: segstart

    !! Initialize the empty matrix
    CALL ConstructEmptyMatrix(ntpolymat, csrmat%num_rows)

    !! We will build the NTPoly matrix from the triplet list.
    CALL ConstructTripletList(triplet_list)

    !! Loop over columns of the CSR matrix that match the NTPoly columns
    iptr = 1
    DO II = ntpolymat%start_column, ntpolymat%end_column-1
       trip%index_column = II
       !! Loop over segments
       DO JJ = csrmat%outer_index(II), csrmat%outer_index(II+1)-1
          segstart = csrmat%segment_index(JJ)
          DO KK = 1, csrmat%segment_len(JJ)
             trip%index_row = segstart + KK - 1
             IF (trip%index_row .GE. ntpolymat%start_row .AND. trip%index_row &
                  & .LT. ntpolymat%end_row) THEN
                trip%point_value = csrmat%values(iptr)
                CALL AppendToTripletList(triplet_list, trip)
             END IF
             iptr = iptr + 1
          END DO
       END DO
    END DO

    !! Preduplicated is set to true because we had each process slice
    !! extract the same information from the CSR mat.
    CALL FillMatrixFromTripletList(ntpolymat, triplet_list, &
         & preduplicated_in=.TRUE.)

    !! Cleanup
    CALL DestructTripletList(triplet_list)

  END SUBROUTINE CSRToNTPoly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert an NTPoly matrix to CSRMat_t.
  !! Note that here I am assuming that the CSRMat is distributed along columns
  !! so you pass in the starting and ending column column of each process.
  SUBROUTINE NTPolyToCSR(ntpolymat, csrmat, start_col, end_col)
    !> The matrix to convert.
    TYPE(Matrix_ps), INTENT(IN) :: ntpolymat
    !> The converted matrix, distributed along columns.
    TYPE(CSRMat_t), INTENT(INOUT) :: csrmat
    !> The starting column
    INTEGER, INTENT(IN) :: start_col
    !> The last column held by a given process
    INTEGER, INTENT(IN) :: end_col
    !! Local Data
    TYPE(TripletList_r) :: triplet_list
    TYPE(Triplet_r) :: trip

    !! First get the triplets associated with the columns of interest
    CALL GetMatrixBlock(ntpolymat, triplet_list, start_row=1, &
         & end_row=csrmat%num_rows, start_column=start_col, &
         & end_column=end_col)

    !! Now convert the triplet list to a csr matrix.
  END SUBROUTINE NTPolyToCSR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE ConversionModule
