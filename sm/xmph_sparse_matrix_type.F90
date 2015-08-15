! Warning: XMPH_GENFILE_COMMENT
!> MAPHYS sparse matrix derived type 
!!
!!----
!! @author Yohan Lee-tin-yien
!! @version v 0.2a 
!!
!! @Date History
!!   Date    : version : Description 
!! - 21/06/11: v0.2b   : Extract MPH_sparse_matrix_enum (Y.Lee-tin-yien)
!! - 13/01/11: v0.2a   : Append comments (Yohan Lee-tin-yien)
!!                       - specifies symmetric cases
!!                       - Append other comments
!! - 2010    : v0.2a   : first implementation (Yohan Lee-tin-yien)
#include "mph_defs_f.h"

Module XMPH_sparse_matrix_type

  ! Dependencies
  Use MPH_sparse_matrix_enum
  Implicit None

  ! 
  Type XMPH_sparse_matrix_t; Sequence

     !> storage format  (see FMT_*      )
     Integer                             :: fmt      
     !> symmetry        (see MATRIX_IS_*)
     Integer                             :: sym      

     !> number of rows
     MPH_INT                          :: m        
     !> number of columns
     MPH_INT                          :: n        
     !> number of non zero
     MPH_INT                          :: nnz      

     !> rows               [nnz]   
     MPH_INT, pointer                 :: i   (:)  
     !> columns            [nnz]   
     MPH_INT, pointer                 :: j   (:)  
     !> values             [nnz]  
     XMPH_FLOAT, pointer               :: v   (:)  

     !> compressed vector  [(ni|nj)+1]
     ! integer, pointer                  :: cs  (:)  
     !> row    compressed vector
     MPH_INT  , dimension(:), pointer :: csr  
     !> column compressed vector
     MPH_INT  , dimension(:), pointer :: csc  

  End Type XMPH_sparse_matrix_t


End Module XMPH_sparse_matrix_type


