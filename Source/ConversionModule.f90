!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This modules has the routines to convert a CSR matrix to NTPoly and
!! back again.
MODULE ConversionModule
  USE FakeChessModule, ONLY : SegMat_t
  USE PSMatrixModule, ONLY : Matrix_ps, GetMatrixTripletList, &
       & ConstructEmptyMatrix, FillMatrixFromTripletList, GetMatrixBlock
  USE SMatrixModule, ONLY : Matrix_lsr
  USE TripletListModule, ONLY : TripletList_r, ConstructTripletList, &
       & AppendToTripletList, DestructTripletList
  USE TripletModule, ONLY : Triplet_r
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: CSRToNTPoly
  PUBLIC :: NTPolyToCSR
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert a SegMat_t type to NTPoly type
  SUBROUTINE CSRToNTPoly(segmat, ntpolymat)
    !> The matrix to convert.
    TYPE(SegMat_t), INTENT(IN) :: segmat
    !> The output matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: ntpolymat
    !! Local Data
    TYPE(segmat_t) :: csrsplit
    TYPE(TripletList_r) :: triplet_list
    TYPE(Triplet_r) :: trip
    INTEGER :: II, JJ, KK
    INTEGER :: iptr
    INTEGER :: segstart

    !! Initialize the empty matrix
    CALL ConstructEmptyMatrix(ntpolymat, segmat%num_rows)

    !! We will build the NTPoly matrix from the triplet list.
    CALL ConstructTripletList(triplet_list)

    !! Loop over columns of the CSR matrix that match the NTPoly columns
    iptr = 1
    DO II = ntpolymat%start_column, ntpolymat%end_column-1
       trip%index_column = II
       !! Loop over segments
       DO JJ = segmat%outer_index(II), segmat%outer_index(II+1)-1
          segstart = segmat%segment_index(JJ)
          DO KK = 1, segmat%segment_len(JJ)
             trip%index_row = segstart + KK - 1
             IF (trip%index_row .GE. ntpolymat%start_row .AND. trip%index_row &
                  & .LT. ntpolymat%end_row) THEN
                trip%point_value = segmat%values(iptr)
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
  !> Convert an NTPoly matrix to segmat_t.
  !! Note that here I am assuming that the segmat is distributed along columns
  !! so you pass in the starting and ending column column of each process.
  SUBROUTINE NTPolyToCSR(ntpolymat, segmat, start_col, end_col)
    !> The matrix to convert.
    TYPE(Matrix_ps), INTENT(IN) :: ntpolymat
    !> The converted matrix, distributed along columns.
    TYPE(SegMat_t), INTENT(INOUT) :: segmat
    !> The starting column
    INTEGER, INTENT(IN) :: start_col
    !> The last column held by a given process
    INTEGER, INTENT(IN) :: end_col
    !! Local Data
    TYPE(Matrix_lsr) :: local_mat
    TYPE(TripletList_r) :: triplet_list, sorted_triplet_list
    TYPE(Triplet_r) :: trip
    INTEGER :: II

    !! First get the triplets associated with the columns of interest
    CALL GetMatrixBlock(ntpolymat, triplet_list, start_row=1, &
         & end_row=segmat%num_rows, start_column=start_col, &
         & end_column=end_col)

    !! Shift the triplets so that they match the data distribution
    DO II = 1, triplet_list%CurrentSize
       triplet_list%data(II)%index_column = &
            & triplet_list%data(II)%index_column - start_col + 1
    END DO

    !! We can now convert the triplet list to csr matrix.

    !! Now convert the triplet list to a csr matrix.
  END SUBROUTINE NTPolyToCSR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE ConversionModule
