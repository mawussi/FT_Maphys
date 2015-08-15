! maphys memory module.
#include "mph_defs_f.h"
#include "mph_macros_f.h"

! [+] module : MPH_mem_mod ---------------------------------------------------
!
!> MaPHyS memory module
!! 
Module MPH_mem_mod
  
  !* No implicit typing *!
  Implicit None

  !* Private constants *!
  Character(len=MAPHYS_STRL), Private, Parameter ::&
       FLNAME = "MPH_mem_mod.F90"

  !* Type definition *! 
  Type MPH_mem_t; sequence

     ! memory usage related to the instance (in Bytes)
     Integer(Kind=8) :: memusage 

     ! memory peak  related to the instance (in Bytes)
     Integer(Kind=8) :: mempeak

     ! memory usage of a package related to the instance (in Bytes).
     ! It is process dependent. 
     Integer(Kind=8) :: pidmemusage 

     ! memory peak of a package related to the instance (in Bytes).
     ! It is process dependent. 
     Integer(Kind=8) :: pidmempeak

  End type MPH_mem_t

  !* Access specifiers *!
  Public :: byte2Mbyte

  Public :: MPH_mem_init
  Public :: MPH_mem_exit
  Public :: MPH_mem_build

  Public :: MPH_mem_setusage
  Public :: MPH_mem_setpeak
  Public :: MPH_mem_setpidusage
  Public :: MPH_mem_setpidpeak

  Public :: MPH_mem_getpeak
  Public :: MPH_mem_getusage
  Public :: MPH_mem_getpidpeak
  Public :: MPH_mem_getpidusage
  Public :: MPH_mem_getnodepeak
  Public :: MPH_mem_getnodeusage

  Public :: MPH_mem_update
  Public :: MPH_mem_add2mem
  Public :: MPH_mem_add2usage
  Public :: MPH_mem_add2peak
  Public :: MPH_mem_add2pidpeak
  Public :: MPH_mem_sub2usage
  
  !* routines *!

  Contains


    ! [+] function : Byte2MByte ------------------------------------------------
    ! 
    !> convert a size in bytes into a size in Mega bytes
    Integer(kind=4) Function Byte2MByte(bytes)
      Implicit None
      Integer(kind=8) :: bytes
      Real   (kind=8) :: conv
      
      conv   = 1_8/Real(1048576,KIND=8)
      Byte2MByte = Ceiling( Real(bytes,KIND=8) * conv )
    End Function Byte2MByte
    
    !
    ! constructors & destructors
    !

  Subroutine MPH_mem_init( mem )
    Type(MPH_mem_t), Intent(out) :: mem
    mem%memusage     = 0
    mem%mempeak      = 0
    mem%pidmempeak   = 0
    mem%pidmemusage  = 0
  End Subroutine MPH_mem_init

  Subroutine MPH_mem_exit( mem )
    Type(MPH_mem_t), Intent(out) :: mem
    mem%memusage     = 0
    mem%mempeak      = 0
    mem%pidmempeak   = 0
    mem%pidmemusage  = 0
  End Subroutine MPH_mem_exit

  Subroutine MPH_mem_build( mem, memusage, mempeak, pidmemusage, pidmempeak )
    Type(MPH_mem_t), Intent(out) :: mem
    Integer(Kind=8) , Intent(in) :: memusage 
    Integer(Kind=8) , Intent(in) :: mempeak
    Integer(Kind=8) , Intent(in) :: pidmemusage 
    Integer(Kind=8) , Intent(in) :: pidmempeak
    mem%memusage     = memusage
    mem%mempeak      = mempeak
    mem%pidmempeak   = pidmemusage
    mem%pidmemusage  = pidmempeak
  End Subroutine MPH_mem_build

  !
  ! setters
  !

  Subroutine MPH_mem_setusage( mem, memset )
    Type(MPH_mem_t), Intent(inout) :: mem
    Integer(Kind=8) , Intent(in) :: memset
    mem%memusage  = memset
    Call MPH_mem_update(mem)
  End Subroutine MPH_mem_setusage

  Subroutine MPH_mem_setpeak( mem, memset )
    Type(MPH_mem_t), Intent(inout) :: mem
    Integer(Kind=8) , Intent(in) :: memset
    mem%mempeak    = memset
    Call MPH_mem_update(mem)
  End Subroutine MPH_mem_setpeak

  Subroutine MPH_mem_setpidusage( mem, memset )
    Type(MPH_mem_t), Intent(inout) :: mem
    Integer(Kind=8) , Intent(in) :: memset
    mem%pidmemusage  = memset
    Call MPH_mem_update(mem)
  End Subroutine MPH_mem_setpidusage

  Subroutine MPH_mem_setpidpeak( mem, memset )
    Type(MPH_mem_t), Intent(inout) :: mem
    Integer(Kind=8) , Intent(in) :: memset
    mem%pidmempeak    = memset
    Call MPH_mem_update(mem)
  End Subroutine MPH_mem_setpidpeak

  !
  ! getters
  !

  Integer(Kind=8) Function MPH_mem_getpidusage(mem)
    Type(MPH_mem_t), Intent(in) :: mem
    MPH_mem_getpidusage = mem%pidmemusage
  End Function  MPH_mem_getpidusage

  Integer(Kind=8) Function MPH_mem_getpidpeak(mem)
    Type(MPH_mem_t), Intent(in) :: mem
    MPH_mem_getpidpeak = mem%pidmempeak
  End Function  MPH_mem_getpidpeak

  Integer(Kind=8) Function MPH_mem_getusage(mem)
    Type(MPH_mem_t), Intent(in) :: mem
    MPH_mem_getusage = mem%memusage
  End Function  MPH_mem_getusage

  Integer(Kind=8) Function MPH_mem_getpeak(mem)
    Type(MPH_mem_t), Intent(in) :: mem
    MPH_mem_getpeak = mem%mempeak
  End Function  MPH_mem_getpeak

  Integer(Kind=8) Function MPH_mem_getallpeak(mem)
    Type(MPH_mem_t), Intent(in) :: mem
    MPH_mem_getallpeak = mem%mempeak + mem%pidmempeak
  End Function  MPH_mem_getallpeak

  Integer(Kind=8) Function MPH_mem_getallusage(mem)
    Type(MPH_mem_t), Intent(in) :: mem
    MPH_mem_getallusage = mem%memusage + mem%pidmemusage
  End Function  MPH_mem_getallusage


  Integer(Kind=8) Function MPH_mem_getnodeusage(mem, env)

    Use mph_error_mod
    Use mph_env_mod
    Implicit None
    Include "mpif.h"

    ! Return the memory usage of the node

    Type(MPH_mem_t), Intent(in) :: mem
    Type(mph_env_t), Intent(in) :: env
    
    Integer :: nodcomm, ierr
    Integer(Kind=8) :: pidmem, nodmem

    nodcomm = mph_env_getNodeComm(env)
    pidmem=mem%memusage + mem%pidmemusage
    Call MPI_AllReduce(pidmem,nodmem,1,&
         MPI_INTEGER8,MPI_SUM,&
         nodcomm, ierr)
    MPH_mem_getnodeusage = nodmem

  End Function  MPH_mem_getnodeusage

  Integer(Kind=8) Function MPH_mem_getnodepeak(mem, env)

    Use mph_env_mod
    Use mph_error_mod
    Implicit None
    Include "mpif.h"

    Type(MPH_mem_t), Intent(in) :: mem
    Type(mph_env_t), Intent(in) :: env
    
    Integer :: nodcomm, ierr
    Integer(Kind=8) :: pidmem, nodmem

    ! Return the memory peak of the node

    nodcomm = mph_env_getNodeComm(env)
    pidmem=mem%mempeak + mem%pidmempeak
    Call MPI_AllReduce(pidmem,nodmem,1,&
         MPI_INTEGER8,MPI_SUM,&
         nodcomm, ierr)
    MPH_mem_getnodepeak = nodmem
  End Function  MPH_mem_getnodepeak

  !
  ! methods
  !

  Subroutine MPH_mem_update(mem)
    Type(MPH_mem_t), Intent(inout) :: mem
    mem%pidmempeak = Max(mem%pidmempeak, mem%pidmemusage)
    mem%mempeak = Max(mem%mempeak, mem%memusage)
  End Subroutine MPH_mem_update

  Subroutine MPH_mem_add2mem( mem, memadded )
    Type(MPH_mem_t), Intent(inout) :: mem
    Type(MPH_mem_t), Intent(in   ) :: memadded

    mem%mempeak     = Max &
         (mem%memusage + memadded%mempeak,mem%mempeak)
    mem%pidmempeak  = Max &
         (mem%pidmemusage + memadded%pidmempeak,mem%pidmempeak)

    mem%memusage    = mem%memusage      + memadded%memusage
    mem%pidmemusage = mem%pidmemusage   + memadded%pidmemusage

    Call MPH_mem_update(mem)
  End Subroutine MPH_mem_add2mem

  Subroutine MPH_mem_add2usage( mem, memadded )
    Type(MPH_mem_t), Intent(inout) :: mem
    Integer(Kind=8) , Intent(in)   :: memadded 
    mem%memusage     = mem%memusage + memadded
    Call MPH_mem_update(mem)
  End Subroutine MPH_mem_add2usage

  Subroutine MPH_mem_add2peak( mem, memadded )
    Type(MPH_mem_t), Intent(inout) :: mem
    Integer(Kind=8) , Intent(in) :: memadded
    Call MPH_mem_update(mem)
    mem%mempeak     = mem%mempeak + memadded
  End Subroutine MPH_mem_add2peak

  Subroutine MPH_mem_add2pidpeak( mem, memadded )
    Type(MPH_mem_t), Intent(inout) :: mem
    Integer(Kind=8) , Intent(in) :: memadded
    Call MPH_mem_update(mem)
    mem%pidmempeak     = mem%pidmempeak + memadded
  End Subroutine MPH_mem_add2pidpeak

  Subroutine MPH_mem_sub2usage( mem, memsubstracted )
    Type(MPH_mem_t), Intent(inout) :: mem
    Integer(Kind=8) , Intent(in)   :: memsubstracted
    mem%memusage     = mem%memusage - memsubstracted
  End Subroutine MPH_mem_sub2usage

  End Module MPH_mem_mod
