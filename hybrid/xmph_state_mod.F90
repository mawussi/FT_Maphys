! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"  
  ! Read user inputs : update ikeep/rkeep, from  current status and icntl/rcntl.
  !
  ! Set of routines to read user inputs.
  ! 
  Module XMPH_state_mod

    Use XMPH_maphys_type
    Use XMPH_state_update_mod
    Use XMPH_state_print_mod
    Use MPH_error_mod
    Implicit None

    ! List of routines
    Public :: XMPH_state_stepstart
    Public :: XMPH_state_stepend

    Private :: XMPH_MAPHYS_Solve_read_options

    ! Private constants
    Character(len=MAPHYS_STRL), Parameter, Private :: &
         FLNAME= "XMPH_ARITHmph_state_mod.F90"


  Contains
    ! [+] routine : XMPH_state_stepstart --------------------------------------
    !
    !> Update the state of the maphys instance at the begining of a step.
    !!
    Subroutine XMPH_state_stepstart & 
         ( mphs, ikeep, rkeep, info )
      Use MPH_maphys_enum
      Use MPH_strat_mod
      Implicit None

      !* Arguments *!

      Type(XMPH_maphys_t), Intent(inout) :: mphs
      Integer       , Intent(inout) :: ikeep  ( MAPHYS_IKEEP_SIZE  ) 
      Real(kind=8)  , Intent(inout) :: rkeep  ( MAPHYS_RKEEP_SIZE  ) 
      Integer       , Intent(  out) :: info


      !- End of header----------------------------------------------------------


      ! set IKEEP_SYMMETRY if IKEEP_CURJOB > CURJOB_Is_Analysis
      If (ikeep(IKEEP_CURJOB) > CURJOB_IsAnalysis) &
           ikeep(IKEEP_SYMMETRY) = mphs%sm_Aii%sym

      ! set the rest
      Call MPH_strat_set(ikeep,rkeep,mphs%icntl,mphs%rcntl,mphs%comm,info)
      MPH_ONFAILURE_ABORT(info)

      ! Inform user
      Call XMPH_state_print_about(mphs)
      Call XMPH_state_print_stepstart(mphs)
      Call XMPH_state_print_listCntl(mphs,info)
      MPH_ONFAILURE_RETURN(info)

      ! Step specific (to refactorize)

      ! Call XMPH_MAPHYS_Solve_read_options(mphs,ikeep,rkeep,info)
      ! MPH_ONFAILURE_RETURN(info)

    End Subroutine XMPH_state_stepstart


    ! [+] routine : XMPH_state_stepend --------------------------------------
    !
    !> Update the state of the maphys instance at the end of a step.
    !!
    Subroutine XMPH_state_stepend(mphs)
      Implicit None
      Type(XMPH_maphys_t), Intent(inout) :: mphs

      Integer :: info_ignored
      
      Call XMPH_state_print_stepend  (mphs)
      Call XMPH_state_UpdateStatus(mphs)
      Call XMPH_state_print(mphs, info_ignored)
      Call XMPH_state_print_listInfo(mphs,info_ignored)
      Call XMPH_state_print_Legend(mphs)

    End Subroutine XMPH_state_stepend


    ! [+] routine : XMPH_MAPHYS_Solve_read_options -------------------------------
    !
    !> Setup the strategy for the solving step according to the context.
    !!
    !! Basically setup the ikeeps and rkeeps for the solving  step,
    !! and check their consistancy.
    !!
    !! @param[in,out] mphs          the maphys instance 
    !! @param[in,out] ikeep         internal controls represented with intege
    !! @param[in,out] rkeep         internal controls represented with reals
    !! @param[   out] info          the routine status
    !!
    !! @author Yohan Lee-tin-yien
    !!
    Subroutine XMPH_MAPHYS_Solve_read_options ( & ! intents
         mphs, ikeep, rkeep,      & ! inout
         info                     & ! out
         )
      Use XMPH_maphys_type
      Use MPH_maphys_enum
      Use mph_log_mod
      Use XMPH_sds_mod, Only : XMPH_sds_t
      Implicit None
      Include "mpif.h"

      !* Subroutine arguments *!

      Type(XMPH_maphys_t), Intent(inout) :: mphs
      Integer       , Intent(inout) :: ikeep  ( MAPHYS_IKEEP_SIZE  ) 
      Real(kind=8)  , Intent(inout) :: rkeep  ( MAPHYS_RKEEP_SIZE  ) 
      Integer       , Intent(  out) :: info

    End Subroutine XMPH_MAPHYS_Solve_read_options




  End Module XMPH_state_mod





