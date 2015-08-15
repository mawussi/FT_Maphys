! Warning: XMPH_GENFILE_COMMENT
Module XMPH_dls_mod

  !* Modules *! 

  Use XMPH_dense_matrix_mod        , Only : &
       XMPH_dense_matrix_t         ! structure
  Use XMPH_dds_mod , Only : &
       XMPH_dds_t  ! structure

  !* No Impilcit typing *! 

  Implicit none

  !* Structure definitions *! 
  
  ! [+] type : dense_linear_system_t -------------------------------------------
  !
  !> We want to solve the dense linear system 
  !! \f[  A . sol = rhs \f]
  !! where A is a dense matrix
  !!
  !! @note About the right-hand-side and solution.
  !! We do not perform in-place computation, even if it is possible in most cases.
  !! The overhead should be acceptable.
  !! This is done especially for Subroutines in mix precision like "DSGESV"
  !!
  Type XMPH_dls_t; Sequence

     !> the dense matrix 
     !!
     !! the dense matrix "A" of the linear system 
     !! - on input, holds  the dense matrix A
     !! - on output, holds the factorisation of A
     !! .
     Type(XMPH_dense_matrix_t)          :: dm_A

     !> the right-hand-side
     !!
     !! Specifies the right-hand-side "rhs" of the linear system.x
     !! - on input, holds the right-hand-side(s)
     !! - on output, may holds the solution
     !! .
     Type(XMPH_dense_matrix_t)          :: dm_B

     !> the solution
     !!
     !! Specifies the solution "sol" of the linear system.x
     !! - on output, holds the solution(s)
     !! .
     Type(XMPH_dense_matrix_t)          :: dm_X

     !> the dense direct solver
     !!
     !! Opaque structure holding the dense direct solver used to solve
     !! the system. (ex : LAPACK)
     !!
     Type(XMPH_dds_t)   :: dds 

     !> the pivots
     !!
     !! Array containing the pivots made during the
     !! resolution of the linear system.
     !! 
     Integer, pointer              :: IPIV (:)
  End Type XMPH_dls_t

End Module XMPH_dls_mod
