! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
! [+] module : XMPH_dds_mod -----------------------------------------------------
!
!> Public interface for a generic dense direct solver
Module XMPH_dds_mod

  !* Module(s) *!
  Use MPH_dds_enum
  Use XMPH_dense_matrix_mod, Only : &
      XMPH_dense_matrix_t        ! derived type

#ifdef HAVE_LIBLAPACK     
  Use XMPH_dds_lapack_mod
#endif

  !* No implicit typing *!
  Implicit none

  !* Derived types *!

  !> Generic object for dense direct solver
  Type XMPH_dds_t; Sequence

     Integer :: choice
#ifdef HAVE_LIBLAPACK     
     Type(XMPH_lapack_t) :: lapack
#endif
  End Type XMPH_dds_t

  !* List of routines *!

  Public :: XMPH_dds_init
  Public :: XMPH_dds_exit
  Public :: XMPH_dds_get_pivots
  Public :: XMPH_dds_factorize
  Public :: XMPH_dds_Solve
  Public :: XMPH_dds_rinfo

Contains

  ! [+] routine :  XMPH_dds_init ------------------------------------------
  !
  !> initialize the instance
  !!
  !! Intialize an instance.
  !!
  !! @param[in,out] dds
  !!
  !!        instance to be initialized
  !!
  !! @param[   out] info
  !!
  !!        routine status
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_dds_init     (dds , info)
    Implicit None

    Type(XMPH_dds_t), intent(inout) :: dds
    integer                    , intent(  out) :: info
    !- End of header -------------------------------------------------

#ifdef HAVE_LIBLAPACK     
    Call XMPH_lapack_init(dds%lapack , info)
#endif

  End Subroutine XMPH_dds_init

  ! [+] routine :  XMPH_dds_exit ------------------------------------------
  !
  !> destroy the dense direct solver
  !!
  !! Deallocate the temporary memory used 
  !! by an instance.
  !!
  !! @param[in,out] dds
  !!
  !!        instance to be freed
  !!
  !! @param[   out] info
  !!
  !!        routine status
  !!
  !! @author Yohan Lee-tin-yien
  !!
  subroutine XMPH_dds_exit (dds, info)
    Implicit None

    ! imported routines :
    Type(XMPH_dds_t), intent(inout) :: dds
    integer                    , intent(  out) :: info
    !- End of header -------------------------------------------------

#ifdef HAVE_LIBLAPACK     
    Call XMPH_lapack_exit (dds%lapack, info)
#endif

  End Subroutine XMPH_dds_exit

  ! [+] routine :  XMPH_dds_get_pivots ----------------------------------------------
  !> return a pointer to the pivots computed during the factorization
  !! 
  !! During factorisation, pivoting may be performed.
  !! This routine return a pointer to those pivots.
  !!
  !! @param[in,out] dds
  !!
  !!        instance of the dense direct solver.
  !!
  !! @param[in,out] piv
  !!
  !!        - On output, points to the pivots.
  !!
  !! @param[   out] info
  !!
  !!        routine status
  !!
  !!
  subroutine XMPH_dds_get_pivots (dds, piv, info)
    Implicit None

    Type(XMPH_dds_t), intent(inout) :: dds
    Integer, pointer           , intent(inout) :: piv (:)
    Integer                    , intent(  out) :: info
    !- End of header -------------------------------------------------
#ifdef HAVE_LIBLAPACK     
    Call XMPH_lapack_get_pivots (dds%lapack, piv,  info)
#endif

  End Subroutine XMPH_dds_get_pivots

  ! [+] routine :  XMPH_dds_factorize -----------------------------------------------
  !> factorize a dense matrix
  !!
  !! Compute the factors of the matrix "MAT".
  !!
  !! During factorisation, pivoting may be performed.
  !! To obtain the pivots, call routine XMPH_dds_get_pivots() .
  !!
  !! @param[in,out] dds
  !!
  !!        instance of the dense direct solver.
  !!
  !! @param[in,out] MAT
  !!
  !!        - On input , holds the matrix to be factorized.
  !!          Its symmetry (MAT%sym) pilots the factorisation type :
  !!          - DM_SYM_isGeneral   -> LU   Factorization
  !!          - DM_SYM_isSPD       -> LLT  Factorization
  !!          - DM_SYM_isSymmetric -> LDLT Factorization
  !!        - On output, holds the factors of the matrix
  !!
  !! @param[   out] info
  !!
  !!        routine status
  !!
  !! @author Yohan Lee-tin-yien
  subroutine XMPH_dds_factorize(dds, MAT, info)
    Implicit None
    Type(XMPH_dds_t), intent(inout) :: dds
    Type(XMPH_dense_matrix_t)       , Intent(inout) :: MAT
    integer                    , intent(  out) :: info
    !- End of header -------------------------------------------------
#ifdef HAVE_LIBLAPACK     
    Call XMPH_lapack_factorize(dds%lapack, MAT, info)
#endif
  End Subroutine XMPH_dds_factorize

  ! [+] routine :  XMPH_dds_solve ---------------------------------------------------
  !> solve the dense linear system
  !!
  !! Solve a dense linear system, knowing 
  !! the factors of the matrix, and the right-hand-side.
  !!
  !! @param[in,out] dds
  !!
  !!        instance of the dense direct solver.
  !!
  !! @param[in    ] MAT
  !!
  !!        holds the factors of the matrix.
  !!        It must be the same MAT of "XMPH_dds_Factorize"
  !!
  !! @param[in    ] RHS
  !!
  !!        - On input , holds the right-hand-side of the linear system.
  !!        - On output, holds the solution        of the linear system.
  !!
  !! @param[   out] info
  !!
  !!        routine status
  !!
  !! @author Yohan Lee-tin-yien
  !!
  subroutine XMPH_dds_solve(dds, MAT, rhs, info)
    Implicit None

    Type(XMPH_dds_t), Intent(inout) :: dds
    Type(XMPH_dense_matrix_t)       , Intent(inout) :: MAT
    Type(XMPH_dense_matrix_t)       , Intent(inout) :: RHS
    integer                    , Intent(  out) :: info
    !- End of header -------------------------------------------------
#ifdef HAVE_LIBLAPACK     
    Call XMPH_lapack_solve(dds%lapack, MAT, rhs, info)
#endif
  End Subroutine XMPH_dds_solve


  ! [+] routine :  XMPH_dds_rinfo ---------------------------------------------------
  !> Obtain an information from a XMPH_dds instance,
  !! stored in a real.
  !!
  !! @param[in    ] dds
  !!
  !!        instance of the dense direct solver.
  !!
  !! @param[in    ] rinfo_flag
  !!
  !!        Flag specifying the wanted informations.
  !!        Possible : all DDS_RINFO*
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Real(kind=8) Function XMPH_dds_rinfo(dds, rinfo_flag )

    !* Arguments *!

    Type(XMPH_dds_t), Intent(in) :: dds
    Integer                    , Intent(in) :: rinfo_flag

    !* Local variables *!
    Real(kind=8), Parameter :: DEFAULT_RES = -9999.d0

    !- End of header -------------------------------------------------

#ifdef HAVE_LIBLAPACK     
    Select Case (rinfo_flag)
    Case(DDS_RINFO_RCOND1); 
       XMPH_dds_rinfo = dds%lapack%rinfo(LAPACK_RINFO_RCOND1)
    Case(DDS_RINFO_NORM1 ); 
       XMPH_dds_rinfo = dds%lapack%rinfo(LAPACK_RINFO_NORM1)
    Case Default          ; 
       XMPH_dds_rinfo = DEFAULT_RES
    End Select
#endif

   End Function XMPH_dds_rinfo



End module XMPH_dds_mod
