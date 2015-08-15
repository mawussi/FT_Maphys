! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"  
  !> Update the state of maphys instance
  Module XMPH_state_update_mod

    ! List of routines 

    Public :: XMPH_state_updateStatus
    Public :: XMPH_state_updateGbStatus

    Private :: XMPH_state_isummarize
    Private :: XMPH_state_rsummarize

    ! Private constants
    Character(len=MAPHYS_STRL), Parameter, Private :: &
         FLNAME= "XMPH_ARITHmph_maphys_read_mod.F90"
    
    Contains


    ! [+] routine : XMPH_state_UpdateStatus -----------------------------------------
    !
    !> update the status of the maphys instance.
    !!
    !! Namely update the informations (iinfo*/rinfo*) of the instance.
    !!  
    !!----
    !!
    !! @param [in,out] mphs
    !!        the maphys instance, it uses the fields:
    !!        - comm        [in ]
    !!        - ikeep       [in ]
    !!        - iinfo*      [in,out]
    !!        - rinfo*      [in,out]
    !!
    !!----
    !!
    !! @author Yohan Lee-tin-yien
    !!
    Subroutine XMPH_state_UpdateStatus(mphs)

      !* Modules & co. *!

      Use MPH_mem_mod
      Use mph_log_mod
      Use mph_error_mod
      Use MPH_maphys_enum
      Use XMPH_maphys_type

      Implicit None

      !* Subroutine arguments *!
      Type(XMPH_maphys_t), Intent(inout) :: mphs

      !* Local Variables *!

      ! Scalars
      Integer    :: iinfo
      Integer    :: rank, master
      Integer    :: sverb,sunit

      !-------------------------------------------------------------------------
      ! [1] Init
      !-------------------------------------------------------------------------

      rank   = mphs%ikeep(IKEEP_MPIRANK)
      master = mphs%ikeep(IKEEP_HOSTRANK)

      !-------------------------------------------------------------------------
      ! [2] Save the memory statistics
      !-------------------------------------------------------------------------

      mphs%iinfo(IINFO_MAPHYS_MEMUSED) = Byte2MByte &
           (MPH_mem_getallusage(mphs%mem))
      mphs%iinfo(IINFO_MAPHYS_MEMPEAK) = Byte2MByte &
           (MPH_mem_getallpeak (mphs%mem))

      mphs%iinfo(IINFO_LIBSDS_MEMUSED) = Byte2MByte &
           (MPH_mem_getpidusage (mphs%mem))
      mphs%iinfo(IINFO_LIBSDS_MEMPEAK) = Byte2MByte &
           (MPH_mem_getpidpeak (mphs%mem))

      mphs%iinfo(IINFO_NODE_MEMUSED) = Byte2MByte &
           (MPH_mem_getnodeusage(mphs%mem,mphs%env))
      mphs%iinfo(IINFO_NODE_MEMPEAK) = Byte2MByte &
           (MPH_mem_getnodepeak (mphs%mem,mphs%env))

      !-------------------------------------------------------------------------
      ! [3] Synchronize the statistics
      !-------------------------------------------------------------------------

      Call XMPH_state_isummarize(&
           mphs%comm, MAPHYS_IINFO_SIZE, mphs%iinfo,&
           mphs%iinfomin,mphs%iinfomax,&
           mphs%iinfoavg,mphs%iinfosig,&
           iinfo)
      CHCKASSRT(iinfo >= 0, iinfo)
      If ( iinfo < 0 ) Goto 9999

      Call XMPH_state_rsummarize(&
           mphs%comm, MAPHYS_RINFO_SIZE, mphs%rinfo,&
           mphs%rinfomin,mphs%rinfomax,&
           mphs%rinfoavg,mphs%rinfosig,&
           iinfo)
      CHCKASSRT(iinfo >= 0, iinfo)
      If ( iinfo < 0 ) Goto 9999

      ! update the global status
      Call XMPH_state_updateGbstatus(mphs)

      !-------------------------------------------------------------------------
      ! Finish
      !-------------------------------------------------------------------------

9999  Continue
      ! set return code
      If ( iinfo < 0 )Then
         mphs%iinfo (1) = iinfo 
         Call XMPH_State_updateGbStatus(mphs)
      End If

    End Subroutine XMPH_state_UpdateStatus


    ! [+] routine : XMPH_state_UpdateGbStatus ----------------------------------
    !
    !> set the global the status of the instance
    !!
    !! Check the status of the instance by
    !! computing iinfog(1) and if an error occured iinfog(2).
    !!
    !!----
    !!
    !! @param [in,out] mphs
    !!        the maphys instance, it uses the fields:
    !!        - comm        [in ]
    !!        - iinfo(1)    [in ]
    !!        - iinfog(1)   [out]
    !!        - iinfog(2)   [in,out]
    !!
    !!----
    !!
    !! @author Yohan Lee-tin-yien
    !!
    !! @todo update
    Subroutine XMPH_state_updateGbStatus(mphs)

      Use XMPH_maphys_type
      Use MPH_maphys_enum, Only : &
           ICNTL_OUTPUT_ErrUnit      ! constant
      Use mph_error_mod
      Use mph_log_mod
      Implicit None
      Include 'mpif.h'

      !* Subroutine arguments *!

      Type(XMPH_maphys_t), Intent(inout) :: mphs

      !* Local Variables *!

      ! Scalars
      Integer :: iinfo
      Integer :: rank, np
      Integer :: i
      Integer :: status
      Integer :: errunit

      ! String
      Character(len=MAPHYS_STRL) :: msg

      ! Array
      Integer, Pointer :: allstatus(:)

      !- End of header----------------------------------------------

      !-------------------------------------------------------------
      ! [1] Gather the local statii
      !-------------------------------------------------------------

      !
      Nullify(allstatus)
      iinfo   = 0
      status  = mphs%iinfo(1)
      errunit = mphs%icntl(ICNTL_OUTPUT_ErrUnit)
      !
      rank = -1
      np   = -1
      !
      Call MPI_Comm_size(mphs%comm, np  , iinfo)
      ASSRT(iinfo == MPI_SUCCESS )
      !
      Call MPI_Comm_rank(mphs%comm, rank, iinfo)
      ASSRT(iinfo == MPI_SUCCESS )
      !
      Allocate( allstatus(np), STAT= iinfo )
      CHCKASSRT( iinfo == 0, iinfo )
      If (iinfo < 0) Goto 9999
      !
      Call MPI_AllGather( &
           status    , 1, MPI_INTEGER, &
           allstatus , 1, MPI_INTEGER, &
           mphs%comm, iinfo )
      ASSRT(iinfo == MPI_SUCCESS)

      !-------------------------------------------------------------
      ! [2] Compute iinfog(1) et iinfog(2)
      !-------------------------------------------------------------
      ! If iinfog(1) = -1 and iinfog(2) corresponds to 
      ! the error on the processor with the smallest rank. 
      mphs%iinfog(1) =  0
      mphs%iinfog(2) =  0
      Do i=1, np
         If ( allstatus(i) < 0 ) Then
            mphs%iinfog(1) = -1
            mphs%iinfog(2) =  i
            Exit
         End If
      End Do

      !-------------------------------------------------------------
      ! [3] Print error message
      !-------------------------------------------------------------

      If ( ( errunit >0) .and. &
           ( rank == 0 ) .and. &
           ( mphs%iinfog(1) == -1 )) Then

         Do i=1, np
            If ( allstatus(i) < 0 ) Then
               Write(msg,*) &
                    "MaPHyS : Error =", allstatus(i), &
                    "on process [", rank+1, "/", np, "]"
               Call mph_log(MSG_ERROR,msg)
            End If
         End Do

      End If

      !-------------------------------------------------------------
      ! [4] Finish
      !-------------------------------------------------------------

9999  Continue
      If (mphs%iinfo(1) < 0 ) Call MPI_Abort(mphs%comm, mphs%iinfo(1), iinfo)
      If (Associated(allstatus)) Deallocate(allstatus)

    End Subroutine XMPH_state_updateGbStatus

  ! [+] routine : XMPH_state_isummarize ----------------------------------------
  !
  !> Get the summarised statistics from an integer statistic.
  !!
  !! @param[in ] comm  The MPI Communicator
  !! @param[in ] nstat The number of statistics
  !! @param[in ] stat  The list of statistic
  !! @param[out] min   The minimum value of the statistic 
  !! @param[out] max   The maximum value of the statistic 
  !! @param[out] avg   The average value of the statistic 
  !! @param[out] sig   The standard deviation of the statistic 
  !! @param[out] info  the routine status
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_state_isummarize &
       ( comm, nstat, istat, imin, imax, iavg, isig, info )
    
    !* Module    *!

    Use mph_error_mod
    Implicit None
    Include "mpif.h"

    !* Arguments *!

    Integer, Intent(in ) ::  comm  
    Integer, Intent(in ) ::  nstat  
    Integer, Intent(in ) ::  istat  (nstat)
    Integer, Intent(out) ::  imin   (nstat)
    Integer, Intent(out) ::  imax   (nstat)
    Real(kind=8), Intent(out) ::  iavg   (nstat)
    Real(kind=8), Intent(out) ::  isig   (nstat)
    Integer, Intent(out) ::  info  

    !* Local variables *!

    ! scalars
    Integer              :: s ! counter of statistics
    Integer              :: np
    Integer              :: i
    Integer              :: val
    Real(kind=8) :: rnp
    Real(kind=8) :: rval

    ! arrays
    Integer, Allocatable :: statlist(:)

    !- End of header -----------------------------------------------------------
    !---------------------------------------------------------------------------
    ! [1] Init 
    !---------------------------------------------------------------------------

    Call MPI_Comm_size( comm, np, info )
    ASSRT( info == MPI_SUCCESS )

    Allocate(statlist(np*nstat),STAT=info )
    CHCKASSRT( info == 0, info )
    If (info < 0) Return

    !---------------------------------------------------------------------------
    ! [2] Gather the values
    !---------------------------------------------------------------------------
    
    Call MPI_AllGather( &
         istat   (1), nstat, MPI_INTEGER, &
         statlist(1), nstat, MPI_INTEGER, &
         comm, info )
    ASSRT( info == MPI_SUCCESS )
    
    rnp = REAL(np,KIND=8)

    !---------------------------------------------------------------------------
    ! [3] Compute min/max/avg 
    !---------------------------------------------------------------------------

    Do s = 1, nstat 

       imin(s) = statlist(s)
       imax(s) = statlist(s)
       iavg(s) = REAL(statlist(s),KIND=8)
       Do i=2,np
          val  = statlist(s+ (i-1)*nstat)
          rval = REAL(val,KIND=8)
          If (imin(s) > val) imin(s) = val
          If (imax(s) < val) imax(s) = val
          iavg(s) = iavg(s) + rval
       End Do
       iavg(s) = iavg(s)/rnp 

    End Do
   
    !---------------------------------------------------------------------------
    ! [4] Compute sig the standard deviation
    !---------------------------------------------------------------------------

    Do s = 1, nstat 

       isig(s) = 0.
       Do i=1,np
          val  = statlist(s+(i-1)*nstat)
          rval = REAL(val,KIND=8)
          isig(s)  = isig(s) + ( rval - iavg(s) )**2
       End Do
       isig(s) = SQRT( isig(s)/rnp )
    End Do

    !---------------------------------------------------------------------------
    ! [5] Finish
    !---------------------------------------------------------------------------

    Deallocate( statlist )

  End Subroutine XMPH_state_isummarize




  ! [+] routine : XMPH_state_rsummarize ----------------------------------------
  !
  !> Get the summarised statistics from an float statistic.
  !!
  !! @param[in ] comm  The MPI Communicator
  !! @param[in ] nstat The number of statistics
  !! @param[in ] rstat  The list of statistics
  !! @param[out] rmin   The minimum value of the statistic 
  !! @param[out] rmax   The maximum value of the statistic 
  !! @param[out] ravg   The average value of the statistic 
  !! @param[out] rsig   The standard deviation of the statistic 
  !! @param[out] info  the routine status
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_state_rsummarize&
       ( comm, nstat, rstat, rmin, rmax, ravg, rsig, info )
    
    !* Module    *!

    Use mph_error_mod
    Implicit None
    Include "mpif.h"

    !* Arguments *!

    Integer, Intent(in ) ::  comm  
    Integer, Intent(in ) ::  nstat 
    Real(kind=8), Intent(in ) ::  rstat (nstat)
    Real(kind=8), Intent(out) ::  rmin  (nstat) 
    Real(kind=8), Intent(out) ::  rmax  (nstat) 
    Real(kind=8), Intent(out) ::  ravg  (nstat) 
    Real(kind=8), Intent(out) ::  rsig  (nstat)
    Integer, Intent(out) ::  info  

    !* Local variables *!

    ! scalars
    Integer              :: np
    Integer              :: p ! counter of processes
    Integer              :: s ! counter on statistics
    Real(kind=8) :: rnp
    Real(kind=8) :: rval

    ! arrays
    Real(kind=8), Allocatable :: statlist(:)

    !- End of header -----------------------------------------------------------
    !---------------------------------------------------------------------------
    ! [1] Init 
    !---------------------------------------------------------------------------

    Call MPI_Comm_size( comm, np, info )
    ASSRT( info == MPI_SUCCESS )

    Allocate(statlist(nstat*np),STAT=info )
    CHCKASSRT( info == 0, info )
    If (info < 0) Return
    
    rnp = REAL(np,KIND=8)

    !---------------------------------------------------------------------------
    ! [2] Gather the values
    !---------------------------------------------------------------------------

    Call MPI_AllGather( &
         rstat(1)   , nstat, MPI_DOUBLE_PRECISION, &
         statlist(1), nstat, MPI_DOUBLE_PRECISION, &
         comm, info )
    ASSRT( info == MPI_SUCCESS )



    !---------------------------------------------------------------------------
    ! [3] Compute min/max/avg 
    !---------------------------------------------------------------------------
    
    Do s=1, nstat

       rmin(s) = statlist(s)
       rmax(s) = statlist(s)
       ravg(s) = statlist(s)
       Do p=2,np
          rval  = statlist(s+(p-1)*nstat)
          If (rmin(s) > rval) rmin(s) = rval
          If (rmax(s) < rval) rmax(s) = rval
          ravg(s) = ravg(s) + rval
       End Do
       ravg(s) = ravg(s)/rnp 
       
    End Do

    !---------------------------------------------------------------------------
    ! [4] Compute sig the standard deviation
    !---------------------------------------------------------------------------
    
    Do s=1, nstat

       rsig(s) = 0.d0
       Do p=1,np
          rval     = statlist(s+(p-1)*nstat)
          rsig(s)  = rsig(s) + ( rval - ravg(s) )**2
       End Do
       rsig(s) = SQRT( rsig(s)/rnp )

    End Do

    !---------------------------------------------------------------------------
    ! [5] Finish
    !---------------------------------------------------------------------------

    Deallocate( statlist )

  End Subroutine XMPH_state_rsummarize


 
  End Module XMPH_state_update_mod
