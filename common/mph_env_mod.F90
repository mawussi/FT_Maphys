#include "mph_defs_f.h"
#include "mph_macros_f.h"

! [+] module : mph_env_mod ---------------------------------------------------
!
!> MaPHyS Environement module
!! 
!! Provides environment detection/manipulation for maphys.
!!
Module mph_env_mod
  
  !* No implicit typing *!
  
  Implicit None

  !* Defined type(s) *!
  
  Include "mph_env_type_f.inc"

  !* Private constants *!
  Character(len=MAPHYS_STRL), Private, Parameter :: FLNAME = &
       "mph_env_mod.F90"

  !* Access specifiers *!

  Public  :: mph_env_init
  Public  :: mph_env_exit
  Public  :: mph_env_setMPIComm
  Public  :: mph_env_getNodeId
  Public  :: mph_env_getNodeComm
  Public  :: mph_env_getNbNodes

  Private :: comm_splitbynode

  !* Routines *!

  Contains
    !
    ! constructors
    ! 

    Subroutine mph_env_init(menv)
      Implicit None
      Type(mph_env_t), Intent(inout) :: menv
      menv%comm     = -1
      menv%nodecomm = -1
      menv%nodeid   = -1
    End Subroutine mph_env_init

    Subroutine mph_env_exit(menv)
      Use mph_error_mod
      Implicit None
      Include "mpif.h"
      Type(mph_env_t), Intent(inout) :: menv

      Logical :: MPICommIsSet 
      Integer :: info
      
      MPICommIsSet = ( menv%nodeid /= -1 )
      If ( MPICommIsSet )Then
         Call MPI_Comm_free(menv%comm,info)
         ASSRT( info == MPI_SUCCESS )

         Call MPI_Comm_free(menv%nodecomm,info)
         ASSRT( info == MPI_SUCCESS )
      End If

      menv%comm     = -1
      menv%nodecomm = -1
      menv%nodeid   = -1
    End Subroutine mph_env_exit
    
    !
    ! setters
    !

    Subroutine mph_env_setMPIComm(menv, comm)

      !* Module(s) *!

      Use mph_error_mod
      Implicit None
      Include "mpif.h"
      
      !* Argument(s) *!

      Type(mph_env_t), Intent(inout) :: menv
      Integer           , Intent(in   ) :: comm

      !* Local variable(s) *!

      Integer :: rank 
      Integer :: info
      Integer :: np

      !- End of header -----------------------------------------------------------

      menv%comm = comm

      Call MPI_Comm_size(menv%comm,np  , info)
      ASSRT( info == MPI_SUCCESS )
      Call MPI_Comm_rank(menv%comm,rank, info)
      ASSRT( info == MPI_SUCCESS )
      Call comm_splitbynode(menv%comm,np,menv%nodecomm,menv%nodeid, info)
      ASSRT( info == MPI_SUCCESS )

#if MPH_ENV_DEBUG
      Call MPI_Comm_size( menv%nodecomm, np, info )
      ASSRT( info == MPI_SUCCESS )

      print *,"NodeId", menv%nodeid, &
           "have ", np, "elements ", &
           "(from Process rank ",rank,")"
#endif


    End Subroutine mph_env_setMPIComm


    !
    ! getters
    !

    Integer Function mph_env_getNodeId(menv)
      Type(mph_env_t), Intent(in) :: menv
      mph_env_getnodeid=menv%nodeid
    End Function mph_env_getnodeid

    Integer Function mph_env_getNodeComm(menv)
      Type(mph_env_t), Intent(in) :: menv
      mph_env_getnodecomm=menv%nodecomm
    End Function mph_env_getnodecomm

    Integer Function mph_env_getNbNodes(menv)
      Use mph_error_mod
      Implicit None
      Include "mpif.h"
      Type(mph_env_t), Intent(in) :: menv
      Integer :: info, res
      Call MPI_AllReduce(&
           menv%nodeid,res,1,MPI_INTEGER,&
           MPI_MAX,menv%comm,info)
      ASSRT( info == MPI_SUCCESS )
      mph_env_getnbnodes = res
    End Function mph_env_getnbnodes

    
    !
    ! Methods
    !
    Subroutine comm_splitbynode(comm, commnp, subcomm, subcommid, info )

      !* Module(s) *!

      Use mph_error_mod
      Implicit None
      Include "mpif.h"
      
      !* Argument(s) *!

      Integer, Intent(in ) :: comm
      Integer, Intent(in ) :: commnp
      Integer, Intent(out) :: subcomm
      Integer, Intent(out) :: subcommid
      Integer, Intent(out) :: info
      
      !* Local variable(s) *!

      ! Scalars
      Integer :: i, st, ed
      Integer :: rank
      Integer :: hkey

      ! Strings
      Character*(MPI_MAX_PROCESSOR_NAME)     :: hname,hname2
      Character*(commnp*MPI_MAX_PROCESSOR_NAME) :: hnames

      !- End of header -----------------------------------------------------------

      hname=""
      hnames=""
      subcomm   = -1
      subcommid = -1
      info      =  MPI_SUCCESS

      Call MPI_Comm_rank(comm,rank, info)
      ASSRT( info == MPI_SUCCESS )

      Call MPI_Get_Processor_name(hname, i ,info)
      ASSRT( info == MPI_SUCCESS )

      Call MPI_AllGather( &
           hname  ,MPI_MAX_PROCESSOR_NAME  ,MPI_CHARACTER, &
           hnames ,MPI_MAX_PROCESSOR_NAME  ,MPI_CHARACTER, &
           comm,info)
      ASSRT( info == MPI_SUCCESS )

      ! get the color
      hname2=hnames(1:MPI_MAX_PROCESSOR_NAME)
      subcommid = 1
      Do i=1,commnp
         st=(i-1)*MPI_MAX_PROCESSOR_NAME + 1
         ed=i *MPI_MAX_PROCESSOR_NAME
         If ( hname2 /= hnames(st:ed)) subcommid = subcommid+1
         If ( hname  == hnames(st:ed)) exit
      End Do

      ! get the key
      hkey=0
      Do i=1,rank+1
         st=(i-1)*MPI_MAX_PROCESSOR_NAME + 1
         ed=i *MPI_MAX_PROCESSOR_NAME
         If ( hnames(st:ed) == hname ) hkey=hkey+1
      End Do

      Call MPI_COMM_SPLIT( comm, subcommid, hkey, subcomm, info)
      ASSRT( info == MPI_SUCCESS )

    End Subroutine comm_splitbynode

End Module mph_env_mod
