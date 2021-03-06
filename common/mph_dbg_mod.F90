#include "mph_defs_f.h"
!> debugging module for maphys
!!
!! It goal is to simply provides
!! shared variables, and/or useful
!! routines for debuging purposes.
!!
!! Namely, provides :
!!
!! - dbg_unit : the unit specifier to write to the debuging file
!! - dbg_rank : the global MPI rank                (ie: from MPI_COMM_WORLD )
!! - dbg_np   : the global MPI number of processes (ie: from MPI_COMM_WORLD )
!! - routines to manipulate the debuging file (modify prefix, etc.)
!!
!! With debugging file in the form of
!!     "dbg_file_prefix","file" 
!! by default
!!     "dbg_file_prefix" = procXXX (where XXX is the dbg_rank+1 )
!!
!! @author Yohan Lee-tin-yien
!!
Module MPH_dbg_mod

  !* No implicit typing *!

  Implicit None

  !* Shared variables *!

  Integer, Public :: dbg_rank
  Integer, Public :: dbg_np
  Integer, Public :: dbg_unit
  Character(MAPHYS_STRL) :: dbg_str
  Character(MAPHYS_STRL) :: dbg_file_prefix

  !* Private variables *!
  Integer, Private, Save :: init_counter = 0

  !* Defined Routines *!

  Public :: MPH_dbg_init
  Public :: MPH_dbg_exit
  Public :: MPH_dbg_set_file
  Public :: MPH_dbg_get_prefix

  Public :: MPH_dbg_prefix_prepend
  Public :: MPH_dbg_prefix_append
  Public :: MPH_dbg_prefix_clear

  Public :: MPH_dbg_barrier

  !* Routines *!

Contains

  !****

  Subroutine MPH_dbg_init()

    Include "mpif.h"

    Integer :: iinfo

    ! Set the shared variables values

    init_counter = init_counter + 1

    If ( init_counter == 1 ) Then
    
       Call MPI_Comm_size(MPI_COMM_WORLD, dbg_np  , iinfo )
       Call MPI_Comm_rank(MPI_COMM_WORLD, dbg_rank, iinfo )

    End If

    Write(dbg_str,'(" proc "I3"/"I3" :")') dbg_rank+1, dbg_np
    dbg_file_prefix = ""
    Write(dbg_file_prefix,'("proc"I3.3".")') dbg_rank+1
    dbg_unit = 11 + dbg_rank

  End Subroutine MPH_dbg_init

  !****

  Subroutine MPH_dbg_exit

    Include "mpif.h"

    Integer :: iinfo

    ! Exit debuging module 
    
    init_counter = init_counter - 1

    Close(UNIT=dbg_unit,IOSTAT=iinfo)

    If ( init_counter <= 0 ) Then
       dbg_unit = -1
       dbg_np   = -1
       dbg_rank = -1
       dbg_str  = ""
       dbg_file_prefix = ""
    End If

  End Subroutine MPH_dbg_exit



  !****

  Subroutine MPH_dbg_set_file(file)

    Character(len=*), intent(in) :: file

    ! Change the dbging output file to "dbg_file_prefix"//"file"

    Close(UNIT=dbg_unit)
    Open (UNIT=dbg_unit, FILE=Trim(dbg_file_prefix)//Trim(file),&
         ACTION="WRITE", STATUS="REPLACE")

  End Subroutine MPH_dbg_set_file

  !****

  Subroutine MPH_dbg_get_prefix (string)

    Character(len=*), Intent(out) :: string

    ! Get the prefix of all debuging files

    string = Trim(dbg_file_prefix)

  End Subroutine MPH_dbg_get_prefix

  !**** 

  Subroutine MPH_dbg_prefix_prepend (string)

    Character(len=*), intent(in) :: string

    ! Preprend "string" to the debuging prefix file 

    dbg_file_prefix = Trim(string) // Trim(dbg_file_prefix)

  End Subroutine MPH_dbg_prefix_prepend

  !****
  
  Subroutine MPH_dbg_prefix_append(string)

    Character(len=*), intent(in) :: string

    ! Append "string" to the debuging prefix file 

    dbg_file_prefix = Trim(dbg_file_prefix) // Trim(string)

  End Subroutine MPH_dbg_prefix_append
  !****

  Subroutine MPH_dbg_prefix_clear()

    ! Delete the debuging prefix file 

    dbg_file_prefix = ""

  End Subroutine MPH_dbg_prefix_clear

  !****

  Subroutine MPH_dbg_barrier()

    include "mpif.h"
    Integer :: ierror

    ! Barrier on MPI_COMM_WORLD
    
    Call MPI_Barrier(MPI_COMM_WORLD, ierror )

  End Subroutine MPH_dbg_barrier

  !****

End Module MPH_dbg_mod
