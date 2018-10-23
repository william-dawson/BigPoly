!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This modules has the routines to convert a CSR matrix to NTPoly and
!! back again.
MODULE ConversionModule
  !! BigDFT Modules
  USE sparsematrix_highlevel, ONLY: ccs_data_from_sparse_matrix
  USE sparsematrix_types, ONLY : matrices, sparse_matrix
  !! NTPoly Modules
  USE sparsematrix_highlevel, ONLY: sparse_matrix_metadata_init_from_file, &
       & sparse_matrix_and_matrices_init_from_file_bigdft, &
       & ccs_data_from_sparse_matrix
  USE PSMatrixModule, ONLY : Matrix_ps, GetMatrixTripletList, &
       & ConstructEmptyMatrix, FillMatrixFromTripletList, GetMatrixBlock
  USE SMatrixModule, ONLY : Matrix_lsr, ConstructMatrixFromTripletList
  USE TripletListModule, ONLY : TripletList_r, ConstructTripletList, &
       & AppendToTripletList, DestructTripletList, SortTripletList
  USE TripletModule, ONLY : Triplet_r
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ChessToNTPoly
  ! PUBLIC :: NTPolyToSeg
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert a CSR matrix type to NTPoly type
  SUBROUTINE ChessToNTPoly(chess_mat, chess_val, ntpolymat)
    TYPE(sparse_matrix):: chess_mat
    TYPE(matrices) :: chess_val
    !> The output matrix.
    TYPE(Matrix_ps), INTENT(INOUT) :: ntpolymat
    !! Local Data
    !! CCS Matrix
    INTEGER, POINTER, DIMENSION(:) :: row_ind, col_ptr
    TYPE(TripletList_r) :: triplet_list
    TYPE(Triplet_r) :: trip
    INTEGER :: II, JJ
    INTEGER :: inner_end
    INTEGER :: iptr
    ! INTEGER :: segstart

    !! Convert to a CCS Matrix
    ALLOCATE(col_ptr(chess_mat%nfvctr))
    ALLOCATE(row_ind(chess_mat%nvctr))
    CALL ccs_data_from_sparse_matrix(chess_mat, row_ind, col_ptr)

    !! Initialize the empty matrix
    CALL ConstructEmptyMatrix(ntpolymat, chess_mat%nfvctr)

    !! We will build the NTPoly matrix from the triplet list.
    CALL ConstructTripletList(triplet_list)

    !! Loop over columns of the CheSS matrix that match the NTPoly columns
    iptr = col_ptr(ntpolymat%start_column)
    DO II = ntpolymat%start_column, ntpolymat%end_column - 1
       !! This early exit condition is needed because of NTPoly's padding.
       IF (II .GT. chess_mat%nfvctr) THEN
          EXIT
       END IF
       !! This check is because chess's CSR datastructure doesn't have the
       !! extract column index showing the total length.
       IF (II .EQ. chess_mat%nfvctr) THEN
          inner_end = chess_mat%nvctr
       ELSE
          inner_end = col_ptr(II+1)-1
       END IF

       !! Loop over values in a row.
       trip%index_column = II
       DO JJ = col_ptr(II), inner_end
          trip%index_row = row_ind(iptr)
          trip%point_value = chess_val%matrix_compr(iptr)
          IF (trip%index_row .GE. ntpolymat%start_row .AND. &
               & trip%index_row .LT. ntpolymat%end_row) THEN
             CALL AppendToTripletList(triplet_list, trip)
          END IF
          iptr = iptr + 1
       END DO
    END DO

    !! Preduplicated is set to true because we had each process slice
    !! extract the same information from the CSR Matrix.
    CALL FillMatrixFromTripletList(ntpolymat, triplet_list, &
         & preduplicated_in=.TRUE.)

    !! Cleanup
    CALL DestructTripletList(triplet_list)

  END SUBROUTINE ChessToNTPoly
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   !> Convert an NTPoly matrix to segmat_t.
  !   !! Note that here I am assuming that the segmat is distributed along columns
  !   !! so you pass in the starting and ending column column of each process.
  !   SUBROUTINE NTPolyToSeg(ntpolymat, segmat, start_col, end_col)
  !     !> The matrix to convert.
  !     TYPE(Matrix_ps), INTENT(IN) :: ntpolymat
  !     !> The converted matrix, distributed along columns.
  !     TYPE(SegMat_t), INTENT(INOUT) :: segmat
  !     !> The starting column
  !     INTEGER, INTENT(IN) :: start_col
  !     !> The last column held by a given process
  !     INTEGER, INTENT(IN) :: end_col
  !     !! Local Data
  !     TYPE(Matrix_lsr) :: csr_mat
  !     TYPE(TripletList_r) :: triplet_list, sorted_triplet_list
  !     TYPE(Triplet_r) :: trip
  !     INTEGER :: II
  !
  !     !! First get the triplets associated with the columns of interest
  !     CALL GetMatrixBlock(ntpolymat, triplet_list, start_row=1, &
  !          & end_row=segmat%num_rows, start_column=start_col, &
  !          & end_column=end_col)
  !
  !     !! Shift the triplets so that they match the data distribution
  !     DO II = 1, triplet_list%CurrentSize
  !        triplet_list%data(II)%index_column = &
  !             & triplet_list%data(II)%index_column - start_col + 1
  !     END DO
  !
  !     !! We can now convert the triplet list to csr matrix.
  !     CALL SortTripletList(triplet_list, segmat%num_columns, &
  !          & segmat%num_rows, sorted_triplet_list)
  !     CALL ConstructMatrixFromTripletList(csr_mat, triplet_list, &
  !          & segmat%num_rows, segmat%num_columns)
  !
  !     !! Now convert the csr_mat to a segmented matrix
  !     CALL CSRToSeg(csr_mat, segmat)
  !   END SUBROUTINE NTPolyToSeg
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   !> Convert a matrix in the compressed sparse row format to a segment
  !   !! matrix.
  !   SUBROUTINE CSRToSeg(csrmat, segmat)
  !     !> Local CSR Matrix to convert.
  !     TYPE(Matrix_lsr) :: csrmat
  !     !> The converted matrix.
  !     TYPE(SegMat_t), INTENT(INOUT) :: segmat
  !     INTEGER II, JJ, iptr
  !
  !     segmat = SegMat_t(csrmat%rows, csrmat%columns)
  !
  !     !! Count the outer index.
  !     iptr = 1
  !     ALLOCATE(segmat%outer_index(segmat%num_columns+1))
  !     segmat%outer_index = 0
  !     segmat%outer_index(1) = 1
  !     DO II = 1, segmat%num_columns
  !        DO JJ = csrmat%outer_index(II), csrmat%outer_index(II+1)
  !           IF (JJ .EQ. 1) THEN
  !              segmat%outer_index(JJ) = segmat%outer_index(JJ) + 1
  !           ELSE IF (csrmat%inner_index(iptr-1) .NE. &
  !                csrmat%inner_index(iptr)) THEN
  !              segmat%outer_index(JJ) = segmat%outer_index(JJ) + 1
  !           END IF
  !           iptr = iptr + 1
  !        END DO
  !     END DO
  !
  !     !! Compute the length of the segments.
  !
  !     !! Copy the data.
  !
  !   END SUBROUTINE CSRToSeg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE ConversionModule
