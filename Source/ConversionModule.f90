!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> This modules has the routines to convert a CSR matrix to NTPoly and
!! back again.
MODULE ConversionModule
  !! BigDFT Modules
  USE sparsematrix_highlevel, ONLY: sparse_matrix_metadata_init_from_file, &
       & sparse_matrix_and_matrices_init_from_file_bigdft, &
       & ccs_data_from_sparse_matrix, &
       & sparse_matrix_init_from_data_ccs, matrices_init_from_data
  USE sparsematrix_types, ONLY : matrices, sparse_matrix
  !! NTPoly Modules
  USE MatrixReduceModule, ONLY : ReduceHelper_t, ReduceAndComposeMatrixData, &
       & ReduceAndComposeMatrixSizes, ReduceAndComposeMatrixCleanup, &
       & TestReduceSizeRequest, TestReduceDataRequest, TestReduceInnerRequest
  USE PSMatrixModule, ONLY : Matrix_ps, GetMatrixTripletList, &
       & ConstructEmptyMatrix, FillMatrixFromTripletList, GetMatrixBlock, &
       & MergeMatrixLocalBlocks
  USE SMatrixModule, ONLY : Matrix_lsr, ConstructMatrixFromTripletList, &
       & DestructMatrix, TransposeMatrix
  USE TripletListModule, ONLY : TripletList_r, ConstructTripletList, &
       & AppendToTripletList, DestructTripletList, SortTripletList
  USE TripletModule, ONLY : Triplet_r
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PUBLIC :: ChessToNTPoly
  PUBLIC :: NTPolyToChess
CONTAINS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert a CheSS matrix type to NTPoly type
  SUBROUTINE ChessToNTPoly(chess_mat, chess_val, ntpolymat)
    !> Input CheSS sparse matrix type
    TYPE(sparse_matrix), INTENT(IN):: chess_mat
    !> Input CheSS matrix of data.
    TYPE(matrices), INTENT(IN) :: chess_val
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Convert an NTPoly matrix to segmat_t.
  !! Note that here I am assuming that the segmat is distributed along columns
  !! so you pass in the starting and ending column column of each process.
  SUBROUTINE NTPolyToChess(ntpolymat, chess_mat, chess_val)
    !> The input matrix to convert.
    TYPE(Matrix_ps), INTENT(INOUT) :: ntpolymat
    !> Output sparse matrix.
    TYPE(sparse_matrix), INTENT(INOUT):: chess_mat
    !> Output matrix of values.
    TYPE(matrices), INTENT(INOUT) :: chess_val
    !! Temporary Variables
    TYPE(Matrix_lsr) :: merged_local_data
    TYPE(Matrix_lsr) :: merged_local_dataT
    TYPE(Matrix_lsr) :: merged_columns
    TYPE(Matrix_lsr) :: merged_columnsT
    TYPE(Matrix_lsr) :: full_gathered
    !! Helpers For Communication
    TYPE(ReduceHelper_t) :: row_helper
    TYPE(ReduceHelper_t) :: column_helper
    !! Helpers for conversion
    INTEGER :: iproc, nproc, comm, nfvctr, nvctr

    !! We gather the NTPoly matrix so it is duplicated on each process.
    !! This is exactly the same logic used in the print matrix routine (for
    !! writing to the console). But I don't wrap it up because there should
    !! be some opportunity for communication hiding here.
    !! Merge all the local data
    CALL MergeMatrixLocalBlocks(ntpolymat, merged_local_data)

    !! Merge Columns
    CALL TransposeMatrix(merged_local_data, merged_local_dataT)
    CALL ReduceAndComposeMatrixSizes(merged_local_dataT, &
         & ntpolymat%process_grid%column_comm, merged_columns, column_helper)
    DO WHILE(.NOT. TestReduceSizeRequest(column_helper))
    END DO
    CALL ReduceAndComposeMatrixData(merged_local_dataT, &
         & ntpolymat%process_grid%column_comm, merged_columns, &
         & column_helper)
    DO WHILE(.NOT. TestReduceInnerRequest(column_helper))
    END DO
    DO WHILE(.NOT. TestReduceDataRequest(column_helper))
    END DO
    CALL ReduceAndComposeMatrixCleanup(merged_local_dataT, merged_columns, &
         & column_helper)

    !! Merge Rows
    CALL TransposeMatrix(merged_columns,merged_columnsT)
    CALL ReduceAndComposeMatrixSizes(merged_columnsT, &
         & ntpolymat%process_grid%row_comm, full_gathered, row_helper)
    DO WHILE(.NOT. TestReduceSizeRequest(row_helper))
    END DO
    CALL ReduceAndComposeMatrixData(merged_columnsT, &
         & ntpolymat%process_grid%row_comm, full_gathered, row_helper)
    DO WHILE(.NOT. TestReduceInnerRequest(row_helper))
    END DO
    DO WHILE(.NOT. TestReduceDataRequest(row_helper))
    END DO
    CALL ReduceAndComposeMatrixCleanup(merged_columnsT, full_gathered, &
          & row_helper)

    !! Convert to CheSS
    iproc = ntpolymat%process_grid%global_rank
    nproc = ntpolymat%process_grid%total_processors
    comm = ntpolymat%process_grid%global_comm
    nfvctr = ntpolymat%actual_matrix_dimension
    nvctr = SIZE(full_gathered%values)
    CALL sparse_matrix_init_from_data_ccs(iproc, nproc, comm, nfvctr, nvctr, &
          & full_gathered%inner_index, full_gathered%outer_index(:nfvctr) + 1, &
          & chess_mat, .FALSE.)
    call matrices_init_from_data(chess_mat, full_gathered%values, chess_val)

    CALL DestructMatrix(merged_local_data)
    CALL DestructMatrix(merged_local_dataT)
    CALL DestructMatrix(merged_columns)
    CALL DestructMatrix(merged_columnsT)
  END SUBROUTINE NTPolyToChess
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE ConversionModule
