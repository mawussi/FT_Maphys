! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
! [+] module : XMPH_sls_mod ----------------------------------------------------
!
!> Module for sparse linear system 
!!
!!
Module XMPH_sls_mod

  Use XMPH_sparse_matrix_type
  Use XMPH_dense_matrix_type
  Use XMPH_sds_mod, Only : XMPH_sds_t

  Implicit None

  !> We want to solve the sparse linear system 
  !! \f[  A . sol = rhs \f]
  !! where A is a sparse matrix
  !!
  !! @par [the sparse matrix]
  !! - sm_A   : the sparse matrix A 
  !!
  !! @par [the right-hand-side]
  !! The right hand side (rhs) can be given in dense or in coordinate format.
  !! - rhs_is_sparse : controls the choosen format of the rhs
  !!          - MPH_TRUE  : the rhs is sparse
  !!          - MPH_FALSE : the rhs is dense
  !!          - Default  = MPH_FALSE
  !!          .
  !! - dm_rhs : pointer to the right-hand-side given in dense format
  !!          It is modified during the solving phase.
  !! - sm_rhs : pointer to the right-hand-side given in coordinate format
  !! . 
  !! @par [the solution]
  !! - dm_sol : a pointer to the solution.
  !!          - with a dense rhs , it points to the dm_rhs 
  !!          - with a sparse rhs, the user initialize/allocate dm_sol
  !!          - before the calling the solving phase
  !!          .
  !! @par [the solver]
  !! - sds    : the sparse direct solver used to solve the system
  !! .
  Type XMPH_sls_t; sequence

     Type(XMPH_sparse_matrix_t)        , pointer :: sm_A
     Type(XMPH_sds_t) , pointer :: sds

     MPH_LOGICAL                         :: rhs_is_sparse
     Type(XMPH_dense_matrix_t)         , pointer :: dm_rhs
     Type(XMPH_sparse_matrix_t)        , pointer :: sm_rhs

     Type(XMPH_dense_matrix_t)         , pointer :: dm_sol

  End Type XMPH_sls_t

  ! list of routines
  Public :: XMPH_SLS_Nullify
  Public :: XMPH_SLS_Free

Contains
  ! [+] routine : XMPH_sls_nullify  --------------------------------------------
  !
  !> nullify a sparse linear system instance
  !! 
  !! @param[in,out ] sls    the sparse linear system instance to nullify
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_SLS_Nullify( sls, info )

    Implicit None

    !* Subroutine arguments *!
    Type(XMPH_sls_t), Intent(inout) :: sls
    Integer, Intent(  out) :: info

    !- End of header -----------------------------------------------------------

    sls%rhs_is_sparse = MPH_LOGICAL_UNSET

    Nullify(sls%sm_A)
    Nullify(sls%sds)
    Nullify(sls%dm_rhs)
    Nullify(sls%sm_rhs)
    Nullify(sls%dm_sol)

    info = 0

  End Subroutine XMPH_SLS_Nullify

  ! [+] routine : XMPH_SLS_Free  -----------------------------------------------
  !
  !> Free a sparse linear system instance
  !! 
  !! Liberate the memory used by a sparse linear system instance.
  !!
  !!-----
  !!
  !! @param[in,out ] sls    the sparse linear system instance to nullify
  !! @param[   out ] info   the routine status
  !!
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  !! @verbatim
  !!  - Date     : Version : Comments
  !!  - 10/01/11 : 0.1a    : Create routine
  !!  - 14/02/11 : 0.1b    : update error/log handling.
  !! @endverbatim
  !!
  Subroutine XMPH_SLS_Free( sls, info )

    !* Modules *!
    Use XMPH_sparse_matrix_mod, Only : &
         XMPH_sm_free ! routine
    Use XMPH_dense_matrix_mod, Only : &
         XMPH_dm_free, & ! routine
         XMPH_dm_Nullify
    Use XMPH_sds_mod, Only : &
         XMPH_sds_Finalize ! routine
    Use mph_log_mod
    Implicit None

    !* Arguments *!
    Type(XMPH_sls_t), Intent(inout) :: sls
    Integer , Intent(  out) :: info

    !* Local variables *!

    ! scalars
    Integer :: iinfo
    Integer :: istep
    Integer :: msgclass

    ! strings
    Character(len=MAPHYS_STRL), Parameter :: rname = "SLS_Free"

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------
    iinfo = 0

    !---------------------------------------------------------------------------
    ! [2] Free the matrix
    !---------------------------------------------------------------------------
    
    If (Associated(sls%sm_A  )) Then
       istep = 2
       Call XMPH_SM_Free(sls%sm_A  , iinfo )
       If (iinfo < 0) Goto 9999
    End If

    !---------------------------------------------------------------------------
    ! [3] Free the solver
    !---------------------------------------------------------------------------

    If (Associated(sls%sds   )) Then
       istep = 3
       Call XMPH_SDS_Finalize (sls%sds  , iinfo )
       If (iinfo < 0) Goto 9999
    End If

    !---------------------------------------------------------------------------
    ! [4] Free the second member 
    !---------------------------------------------------------------------------

    Select Case (sls%rhs_is_sparse)
    Case (MPH_TRUE) 

       If (Associated(sls%sm_rhs )) Then
          istep = 41
          Call XMPH_SM_Free   (sls%sm_rhs, iinfo )
          If (iinfo < 0) Goto 9999
       End If

       If ( Associated(sls%dm_sol) ) Then
          istep = 42
          Call XMPH_DM_Free (sls%dm_sol, iinfo )
          If (iinfo < 0) Goto 9999
       End If

    Case (MPH_FALSE)

       If ( Associated(sls%dm_rhs) ) Then
          istep = 43
          Call XMPH_DM_Free    (sls%dm_rhs, iinfo )
          If (iinfo < 0) Goto 9999
       End If

       If ( Associated(sls%dm_sol) ) Then
          istep = 44
          Call XMPH_DM_Nullify (sls%dm_sol, iinfo )
          If (iinfo < 0) Goto 9999
       End If

    Case Default ! Do Nothing
       Continue
    End Select


    !---------------------------------------------------------------------------
    ! [5] Free the components
    !---------------------------------------------------------------------------

    If ( Associated(sls%sm_A  ) ) Deallocate(sls%sm_A)
    If ( Associated(sls%sds   ) ) Deallocate(sls%sds)
    If ( Associated(sls%dm_rhs) ) Deallocate(sls%dm_rhs)
    If ( Associated(sls%sm_rhs) ) Deallocate(sls%sm_rhs)
    If ( Associated(sls%dm_sol) ) Deallocate(sls%dm_sol)

    !---------------------------------------------------------------------------
    ! [6] Finish
    !---------------------------------------------------------------------------

9999 Continue
  ! Print error/warning messages
  If ( iinfo /=  0 ) Then
     
     If ( iinfo > 0) msgclass = MSG_WARNING
     If ( iinfo < 0) msgclass = MSG_ERROR
     
     Select Case(istep) 
     Case Default
        Call mph_logWithInfo (msgclass,istep, Trim(rname)//&
             " at internal step =")
     End Select
     
  End If

  ! Set return code
  If ( iinfo == 0 ) info =  0
  If ( iinfo <  0 ) info = -istep
  If ( iinfo >  0 ) info = +istep


End Subroutine XMPH_SLS_FREE




End Module XMPH_sls_mod
