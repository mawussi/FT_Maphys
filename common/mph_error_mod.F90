#include "mph_defs_f.h"
!> Module to print maphys messages & handle errors 
!!
!! @note
!! The routines MPH_err_init() and MPH_err_exit() are not thread-safe.
Module mph_error_mod

  !* No Implicit Typing *!
  Implicit None

  !* Private variables *!
  Integer, Private, Save :: errcomm   = -1

  !* Access specifiers to routines *!

  Public :: mph_err_init
  Public :: mph_err_exit
  Public :: mph_abort
  Public :: mph_assert
  Public :: mph_check

  !* Routines *!

Contains

  !> Initialize the module.
  !!
  !! Set the communicator to use for errors.
  !!
  !! @param[in] comm      MPI communicator used for errors (MPI_Abort)
  !! 
  !! @note 
  !! This routine would be called at the begining of maphys_driver.
  !! (ie, each time the user call maphys)
  Subroutine mph_err_init ( comm )

    Integer, Intent(in) ::  comm

    errcomm   = comm

  End Subroutine mph_err_init

  !> Exit the module.
  !!
  !! Reset private variables to their default values.
  !! 
  Subroutine mph_err_exit

    errcomm   = -1

  End Subroutine mph_err_exit

  !> Abort a maphys instance
  Subroutine mph_abort

    !* Module(s) & co. *!

    Implicit None
    Include "mpif.h"   
         
    !* Local variables.
    Integer :: iinfo

    Call MPI_Abort(errcomm, -1, iinfo)

  End Subroutine mph_abort

  
  !> Check the return of MPI.
  !!
  !! @param [in,out] iinfo 
  !!   - On input , holds the return call to an MPI routine
  !!   - On output, holds 0 if MPI routine succeeded, -1, if not.
  !!                     
  Subroutine mph_error_CheckMPI(iinfo)

    !* Module(s) & co. *!

    Implicit None
    Include "mpif.h"   
         
    !* Local variables.
    Integer, Intent(out) :: iinfo

    !- End of header -----------------------------------------------------------

    If ( iinfo /= MPI_SUCCESS ) Then
       iinfo = -1
    Else
       iinfo = 0
    End If

  End Subroutine mph_error_CheckMPI


  !> Assert a statement
  Subroutine mph_assert(assertion, line, flname )

    !* Module(s) *!

    Use mph_log_mod
    Implicit None

    !* Arguments *!
    Logical, Intent(in ) :: assertion
    Integer, Intent(in ) :: line
    Character(len=*), Intent(in) :: flname

    !* Local variables *!
    ! strings
    Character(len=MAPHYS_STRL) :: msg

    !- End of header -----------------------------------------------------------
    
    If ( assertion .eqv. .False. ) Then
       Write(msg,'(2A,I6,A)') Trim(flname),":line",line,". Assertion Failed"
       Call mph_log(MSG_ERROR,msg)
       Call mph_abort
    End If
  
  End Subroutine mph_assert
  
  !> Check that a statement is True 
  !!
  !! if the "statement" is false print useful message.
  !! and set "info" to "-line"
  !! 
  !! @param [in] assertion   the assertion to check
  !! @param [in] line        the line number (should be __LINE__)
  !! @param [inout] info     on success, do nothing,
  !!                         on failure it is set to -line
  !! @param [in] flname      the filename
  Subroutine mph_check(statement, line, info, flname )
    
    !* Module(s) *!

    Use mph_log_mod
    Implicit None

    !* Arguments *!
    Logical, Intent(in ) :: statement
    Integer, Intent(in ) :: line
    Integer, Intent(inout) :: info
    Character(len=*), Intent(in) :: flname

    !* Local variables *!
    ! strings
    Character(len=MAPHYS_STRL) :: msg
    
    !- End of header -----------------------------------------------------------
    
    If ( .Not. statement ) Then
       Write(msg,'(2A,I6,A)') Trim(flname),":line",line,". Check Failed"
       Call mph_log(MSG_ERROR,msg)
       info = - line
    End If

  End Subroutine mph_check

End Module mph_error_mod
