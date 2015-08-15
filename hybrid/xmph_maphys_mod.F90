! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"
!#define MAPHYS_DEBUG 1
#if MAPHYS_DEBUG
#define DBG(x) x
#else
#define DBG(x) ! x
#endif
  
  !> MAPHYS top module
  !! 
  !! It is basically the "main" of the library.
  !! It handles :
  !!   - the calls of the routines according to the resolution strategy
  !!   - the statistics (timers, memory counters and backward errors)
  !!
  Module XMPH_maphys_mod

    !* Module(s) *!
    Use MPH_error_mod
    Use MPH_time_mod
    Use MPH_maphys_enum
    Use XMPH_maphys_type
    Use XMPH_maphys_aux_mod

    !* No implicit typing *!
    Implicit None

    !* Private constants *!
    Character(len=MAPHYS_STRL), Private, Parameter :: &
         FLNAME= "XMPH_ARITHmph_maphys_mod.F90"

    !* Defined routine(s) *!
    Public :: XMPH_MAPHYS_Driver
    Public :: XMPH_MAPHYS_Init
    Public :: XMPH_MAPHYS_Analyze
    Public :: XMPH_MAPHYS_Factorize
    Public :: XMPH_MAPHYS_Precond
    Public :: XMPH_MAPHYS_Solve
    Public :: XMPH_MAPHYS_Exit
    Private :: XMPH_maphys_pretreatInputMatrix
    Private :: XMPH_EstimateError

  Contains

    ! [+] routine : XMPH_MAPHYS_Driver  ----------------------------------------
    !
    !> This is the maphys driver.
    !!
    !! This routine, launch the appropriate step(s) in the algorithm, 
    !! and perform global error handling.
    !!
    !! 
    Subroutine XMPH_maphys_Driver(mphs)

      !* Module(s) & co *!

      Use mph_log_mod
      Use XMPH_state_mod
      Implicit None
      include 'mpif.h'

      !* Arguments *!

      Type(XMPH_maphys_t), intent(inout) :: mphs

      !* Local variables *!

      ! Scalars

      Integer :: iinfo
      Integer :: rank
      Logical :: is_MPI_initialized

      !- End of header--------------------------------------------------------------

      !-----------------------------------------------------------------------------
      ! [-] Initialize local variables
      !-----------------------------------------------------------------------------
      !

      is_MPI_initialized = .False.
      iinfo = MPH_SUCCESS

      !-----------------------------------------------------------------------------
      ! [-] Verify that a few parameters are set
      !-----------------------------------------------------------------------------
      !
      ! MPI 
      mphs%iinfo( IINFO_STATUS ) = 0
      Call MPI_Initialized(is_MPI_initialized, iinfo)
      CHCKASSRT( iinfo == MPI_SUCCESS, mphs%iinfo( IINFO_STATUS ))

      Call MPI_Comm_rank(mphs%comm,rank,iinfo)
      CHCKASSRT( iinfo == MPI_SUCCESS, mphs%iinfo( IINFO_STATUS ))

      !-----------------------------------------------------------------------------
      ! [-] Init several modules
      !-----------------------------------------------------------------------------

      Call mph_err_init( mphs%comm )
      Call mph_log_init(&
           mphs%icntl(ICNTL_OUTPUT_ErrUnit),mphs%icntl(ICNTL_OUTPUT_WarnUnit),&
           mphs%icntl(ICNTL_OUTPUT_StdUnit),mphs%icntl(ICNTL_OUTPUT_Verbosity),&
           rank)

      !-----------------------------------------------------------------------------
      ! [-] Select the step to launch
      !-----------------------------------------------------------------------------
      !
      Select Case(mphs%job)
      Case (-1)
         Call XMPH_MAPHYS_Init (mphs)
         Call XMPH_state_updateGbstatus(mphs)
      Case (-2)
         Call XMPH_MAPHYS_Exit (mphs)
         Call XMPH_state_updateGbstatus(mphs)
      Case (1)
         Call XMPH_MAPHYS_Analyze (mphs)
      Case (2) 
         Call XMPH_MAPHYS_Factorize (mphs)
      Case (3) 
         Call XMPH_MAPHYS_Precond (mphs)
      Case (4) 
         Call XMPH_MAPHYS_Solve (mphs)
      Case (6)
         Call XMPH_MAPHYS_Analyze (mphs)
         If ( mphs%iinfog(IINFOG_STATUS) < 0) Goto 9999

         Call XMPH_MAPHYS_Factorize (mphs)
         If ( mphs%iinfog(IINFOG_STATUS) < 0) Goto 9999

         Call XMPH_MAPHYS_Precond (mphs)
         If ( mphs%iinfog(IINFOG_STATUS) < 0) Goto 9999

         Call XMPH_MAPHYS_Solve (mphs)
         If ( mphs%iinfog(IINFOG_STATUS) < 0) Goto 9999
      Case (7)
         Call XMPH_MAPHYS_Factorize (mphs)
         If ( mphs%iinfog(IINFOG_STATUS) < 0) Goto 9999

         Call XMPH_MAPHYS_Precond (mphs)
         If ( mphs%iinfog(IINFOG_STATUS) < 0) Goto 9999
      End Select

      !-------------------------------------------------------------------------
      ! [-] Exit several modules
      !-------------------------------------------------------------------------
9999  Continue
      Call mph_err_exit
      Call mph_log_exit

    end Subroutine XMPH_maphys_Driver


    ! [+] routine : MAPHYS_Init  -----------------------------------------------
    !
    !> initialize a maphys_instance
    !! 
    !! @param[in,out ] mphs     the maphys instance to initialize
    !!
    !! @author Yohan Lee-tin-yien
    !!
    Subroutine XMPH_maphys_Init(mphs) ! inout

      !* Modules & co *!
      Use mph_log_mod
      Use MPH_mem_mod
      Use mph_env_mod
      Use MPH_strat_mod
      Implicit None
      Include "mpif.h"

      !* Arguments *!

      Type(XMPH_maphys_t), Intent(inout) :: mphs

      !* Local variables *!

      ! Scalars
      Integer                      :: iinfo
      Integer                      :: rank,np

      !- End of header ---------------------------------------------------------

      !-------------------------------------------------------------------------
      ! [-] Initialize local variables
      !-------------------------------------------------------------------------

      iinfo = MPH_SUCCESS
      mphs%ikeep(IKEEP_CURJOB) = CURJOB_IsInit
      
      ! temporaly set the logging error to standard output.

      Call mph_log_init(0,0,6,5,-1)
      Call MPI_Comm_size(mphs%comm, np  , iinfo )
      ASSRT(iinfo == MPI_SUCCESS)
      Call MPI_Comm_rank(mphs%comm, rank, iinfo )
      ASSRT(iinfo == MPI_SUCCESS)
      Call mph_log_init(6,6,6,5,rank)

      !-------------------------------------------------------------------------
      ! [-] Nullify all components & Give default values to parameters
      !-------------------------------------------------------------------------

      Call XMPH_MAPHYS_Nullify(mphs,iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)

      Call XMPH_MAPHYS_Set_default_icntl(mphs%icntl,iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)

      Call XMPH_MAPHYS_Set_default_rcntl(mphs%rcntl,iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)

      !-------------------------------------------------------------------------
      ! [-] Set other components
      !-------------------------------------------------------------------------

      Write(mphs%version,*,IOSTAT=iinfo) MAPHYS_VERSION
      CHCKASSRT(iinfo == 0, iinfo )
      If (iinfo < 0) Goto 9999

      Call MPH_strat_init(mphs%ikeep,mphs%rkeep)
      
      Call MPH_mem_init(mphs%mem)
      Call mph_env_init(mphs%env)
      Call mph_env_setMPIComm(mphs%env,mphs%comm)

      mphs%ikeep(IKEEP_MPIRANK)  = rank
      mphs%ikeep(IKEEP_HOSTRANK) = 0
      mphs%ikeep(IKEEP_NBDOMAINS) = np
      mphs%ikeep(IKEEP_NODES   ) = mph_env_getNbNodes (mphs%env)
      mphs%ikeep(IKEEP_NODEID  ) = mph_env_getNodeId  (mphs%env)
      mphs%ikeep(IKEEP_NODECOMM) = mph_env_getNodeComm(mphs%env)

      mphs%iinfo(IINFO_NODEID)   = mphs%ikeep(IKEEP_NODEID)
      Call MPI_Comm_size(&
           mph_env_getNodeComm(mphs%env), &
           mphs%iinfo(IINFO_NODEDOMAINS),&
           iinfo)
      ASSRT( iinfo == MPI_SUCCESS )

      !-------------------------------------------------------------------------
      ! [*] Finish
      !-------------------------------------------------------------------------
      
9999  Continue
      mphs%iinfo( IINFO_STATUS ) =  iinfo

    end Subroutine XMPH_maphys_Init


    ! [+] routine : MAPHYS_Analyze ---------------------------------------------
    !> perform the analysis phase
    !!
    !!  distribute the global linear system into local linear systems (one
    !!  per MPI process).
    !!
    !!  @param[in,out] mphs the maphys instance
    !!
    !!  @author Yohan Lee-tin-yien
    !!
    Subroutine XMPH_maphys_Analyze (mphs)

      !* Modules & co. *!
      Use MPH_mem_mod
      Use MPH_log_mod
      Use XMPH_distsystem_mod
      Use XMPH_sparse_matrix_mod
      Use XMPH_state_mod      

      Implicit None
      Include 'mpif.h'

      !* Arguments *!
      Type(XMPH_maphys_t), Intent(inout) :: mphs

      !* Local Variables *!

      ! Scalars
      Integer :: iinfo, ignore
      Integer    :: rank, master
      Integer :: verbosity      ! maphys level of verbosity
      Integer, Parameter :: MSGLEN= 2048
      Character(len=MSGLEN) :: msg ! a message.
      ! Derived types
      Type(XMPH_sparse_matrix_t) :: gb_A

      !- End of header----------------------------------------------------------

      !-------------------------------------------------------------------------
      ! [-] Init
      !-------------------------------------------------------------------------

      ! [--] Scalars
      iinfo = MPH_SUCCESS
      CHCKASSRT( mphs%iinfo( IINFO_STATUS ) >= 0, iinfo)
      MPH_ONFAILURE_RETURN(iinfo)

      mphs%ikeep(IKEEP_CURJOB) = CURJOB_IsAnalysis

      ! [--] Timers
      Call MPH_time_start(mphs%rinfo (RINFO_TIMING_Total))
      Call MPH_time_start(mphs%rinfo (RINFO_TIMING_Analysis))

      ! [--] Get MPI information
      rank   = mphs%ikeep(IKEEP_MPIRANK)
      master = mphs%ikeep(IKEEP_HOSTRANK)

      ! [--] Controls
      Call XMPH_state_stepstart(mphs, mphs%ikeep, mphs%rkeep,iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)
      !CHCKASSRT( .False., ignore)

      !-------------------------------------------------------------------------
      ! [-] Get the distributed system
      !-------------------------------------------------------------------------

      ! [--] Get the local matrix, domain, rhs partition

      Select Case (mphs%ikeep(IKEEP_INSYSTEM))
      Case(INSYSTEM_IsCentralized)

      !   CHCKASSRT( .False., ignore)
         Call XMPH_maphys_pretreatInputMatrix(mphs,master,rank,gb_A, iinfo)
         MPH_ONFAILURE_GOTO9999(iinfo)

      !   CHCKASSRT( .False.,ignore)
         Call XMPH_distsystem_compute(mphs,master,gb_A,mphs%comm,iinfo)
         MPH_ONFAILURE_GOTO9999(iinfo)

      !   CHCKASSRT( .False.,ignore)
         Call XMPH_sm_free(gb_A, iinfo)
         MPH_ONFAILURE_GOTO9999(iinfo)

      !   CHCKASSRT( .False.,ignore)
         Call MPH_mem_sub2usage(mphs%mem, XMPH_sm_sizeof( gb_A ))

      Case(INSYSTEM_IsDistributed)
         
      !   CHCKASSRT( .False., ignore)
         Allocate(mphs%sls_domain%sm_A, STAT=iinfo)
         CHCKALLOC(iinfo)
         MPH_ONFAILURE_GOTO9999(iinfo)

         Call mph_log_GetVerbosity(verbosity)
         If ( verbosity >= MSG_VERBOSE ) Then
            Write(msg,'(2A,I5,A,I10)') 'Domain Decomposition ->',&
                 ' local interface [',rank,'] -> size = ',&
                 mphs%lc_domain%myndofintrf
            Call mph_log(MSG_VERBOSE, msg)
            Write(msg,'(2A,I5,A,I10)') 'Domain Decomposition ->',&
                 ' local interior [',rank,'] -> size = ',&
                 mphs%lc_domain%myndofinterior
            Call mph_log(MSG_VERBOSE, msg)
         End If
         
         Call XMPH_sm_CreateFromData(&
              mphs%sls_domain%sm_A,SM_FMT_IJV,mphs%sym,&
              mphs%n, mphs%n, mphs%nnz, &
              mphs%rows, mphs%cols, mphs%values, iinfo )
         MPH_ONFAILURE_GOTO9999(iinfo)

      Case Default

         CHCKASSRT(.False.,iinfo)
         MPH_ONFAILURE_GOTO9999(iinfo)

      End Select
      
      ! [--] Build the subblocs 

      Call XMPH_distsystem_buildSubBlocs(mphs, iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)

      !-------------------------------------------------------------------------
      ! [-] Exit
      !-------------------------------------------------------------------------

9999  Continue

      ! [--] timers

      Call MPH_time_stop(mphs%rinfo (RINFO_TIMING_Analysis))

      ! [--] set return code
      mphs%iinfo( IINFO_STATUS ) = iinfo 
      Call XMPH_state_stepend(mphs)

    End Subroutine XMPH_maphys_Analyze


    ! [+] routine : MAPHYS_Factorize--------------------------------------------
    !
    !> perform the factorization step
    !! 
    !! Compute the factors of each local system + extract the schur complement
    !!
    !! @param[in,out] mphs   the maphys instance
    !!
    !! @author Luc   Giraud
    !! @author Azzam Haidar
    !! @author Yohan Lee-tin-yien
    !!
    Subroutine XMPH_maphys_Factorize(mphs)

      !* Modules & co. *!

      Use mph_log_mod
      Use XMPH_domainsls_mod
      Use XMPH_state_mod
      Use XMPH_dense_matrix_mod, Only :&
           XMPH_dm_sizeof
      Use XMPH_sparse_matrix_mod, Only : &
           XMPH_sm_sizeof
      Use XMPH_sds_mod, Only : & ! routine(s)
           XMPH_sds_get_numstats, &
           XMPH_sds_get_memstats, &
           XMPH_sds_get_perfstats
      Use MPH_mem_mod
      DBG(Use mph_dbg_mod)

      Implicit None
      Include "mpif.h"

      !* Arguments *!

      Type(XMPH_maphys_t), Intent(inout) :: mphs

      !* Local Variables *!
      
      ! Constants
      Integer, Parameter :: UNUSED_INT = -1
      Real(KIND=XMPH_FLOATKIND), Parameter :: UNUSED_FLOAT = XMPH_FLOATZERO

      ! Scalars
      Integer    :: iinfo
      Integer    :: comm
      Integer    :: schur_strategy
      Integer    :: which_sds
      MPH_INT :: nschur 
      Integer(kind=8) :: piv

      MPH_INT :: pilut_LUnzRowMax
      MPH_INT :: pilut_SnzRowMax
      Real(KIND=XMPH_FLOATKIND) :: pilut_LUtol
      Real(KIND=XMPH_FLOATKIND) :: pilut_Stol

      ! Arrays
      Integer :: thread_icntl(MTHREAD_ICNTL_SIZE)
     
      ! String
      Character(len=MAPHYS_STRL) :: msg

      ! Derived types
      Type(XMPH_sparse_matrix_t       ), Pointer :: A
      Type(XMPH_sds_t), Pointer :: Aii_factors
      Type(MPH_mem_t)                    :: sds_mem

      !- End of header----------------------------------------------------------

      !-------------------------------------------------------------------------
      ! [-] Init
      !-------------------------------------------------------------------------

      Nullify(Aii_factors,A)
      iinfo = MPH_SUCCESS
      CHCKASSRT( mphs%iinfo( IINFO_STATUS ) >= 0, iinfo)
      MPH_ONFAILURE_RETURN(iinfo)

      mphs%ikeep(IKEEP_CURJOB) = CURJOB_IsFacto
      Call MPH_time_start(mphs%rinfo(RINFO_TIMING_Facto))

      !-------------------------------------------------------------------------
      ! [--] Select the strategy 
      !-------------------------------------------------------------------------

      Call XMPH_state_stepstart(mphs, mphs%ikeep, mphs%rkeep,iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)

      !-------------------------------------------------------------------------
      ! [--] Initialize scalars, arrays, etc.
      !-------------------------------------------------------------------------

      schur_strategy  = mphs%ikeep(IKEEP_SCHUR_STRATEGY)
      which_sds       = mphs%ikeep(IKEEP_SDS_Facto)

      comm            = MPI_COMM_SELF
      nschur          = mphs%lc_domain%myndofintrf

      pilut_LUnzRowMax  = mphs%ikeep( IKEEP_ILU_LUFill )  
      pilut_SnzRowMax   = mphs%ikeep( IKEEP_ILU_SCHURFill )   
      pilut_LUtol    = mphs%rkeep( RKEEP_ILU_LUThreshold )    
      pilut_Stol     = mphs%rkeep( RKEEP_ILU_SCHURThreshold  )    

      ! Arrays
      thread_icntl(1) = mphs%icntl(ICNTL_2LVLS_Bind)
      thread_icntl(2) = mphs%icntl(ICNTL_2LVLS_NNodes)
      thread_icntl(3) = mphs%icntl(ICNTL_2LVLS_NcoresPerNode)
      thread_icntl(4) = mphs%icntl(ICNTL_2LVLS_NThrdsPerProc)
      thread_icntl(5) = mphs%icntl(ICNTL_2LVLS_NProcs)

      ! Derived types

      Allocate(mphs%sls_domain%sds)
      Aii_factors => mphs%sls_domain%sds
      A => mphs%sls_domain%sm_A

      !-----------------------------------------
      ! [-] Factorize Aii and compute the Schur 
      !-----------------------------------------

      Call MPH_Time_start(mphs%rinfo(RINFO_TIMING_FacLcMat))
      Call MPH_Time_start(mphs%rinfo(RINFO_TIMING_GetSchur))

      Select Case(schur_strategy)

         !----------------------------------------------------------------
         ! [--] Get the schur & the Aii's factors from SDS
         !----------------------------------------------------------------

      Case (SCHUR_STRATEGY_Exact);

         Call XMPH_domainsls_GetdmSchurAndFacto(      & ! intents
              comm, thread_icntl,which_sds, & ! in
              nschur, A ,                   & !
              Aii_factors,                  & ! inout
              mphs%dm_schur,  iinfo         & ! out
              )
         MPH_ONFAILURE_GOTO9999(iinfo)

         DBG(Call MPH_dbg_init)
         DBG(Call MPH_dbg_set_file("dmschurcoo.mtx"))
         DBG(Call XMPH_dm_mmwritecoo(mphs%dm_schur, dbg_unit, iinfo ))
         DBG(Call MPH_dbg_exit)

         !----------------------------------------------------------------
         ! [--] Estimates the schur with PILUT &
         !       the Aii's factors from SDS  
         !----------------------------------------------------------------

      Case (SCHUR_STRATEGY_ApprxWithIluT);

         Call XMPH_domainsls_GetSmSchurAndFacto(    & 
              comm, thread_icntl, which_sds, &
              mphs%sm_Aii, Aii_factors, nschur, A, &
              pilut_LUnzRowMax,pilut_LUtol,&
              pilut_SnzRowMax,pilut_Stol, & 
              mphs%sm_schur,  iinfo )
         MPH_ONFAILURE_GOTO9999(iinfo)

         DBG(Call MPH_dbg_init) 
         DBG(Call MPH_dbg_set_file("smschur.mtx"))
         DBG(Call XMPH_sm_mmwrite(mphs%sm_schur, dbg_unit, iinfo))
         DBG(Call MPH_dbg_exit)

         !----------------------------------------------------------------
         ! [--] Default case
         !----------------------------------------------------------------

      Case Default;

         Write(msg,"(A,I2,A)") &
              "Bad internal value for IKEEP(",IKEEP_SCHUR_STRATEGY,")"
         Call MPH_Log(MSG_ERROR, msg )

         Write(msg,"(A,I2,A)") &
              "Origin : ICNTL(",ICNTL_SCHUR_STRATEGY,") may be wrong" 
         Call MPH_Log(MSG_ERROR, msg )
         
         CHCKASSRT(.False., iinfo)
         MPH_ONFAILURE_GOTO9999(iinfo)

      End Select

      Call MPH_Time_stop(mphs%rinfo(RINFO_TIMING_FacLcMat))
      Call MPH_Time_stop(mphs%rinfo(RINFO_TIMING_GetSchur))
      Call MPH_time_stop(mphs%rinfo(RINFO_TIMING_Facto))

      !-------------------------------------------------------------------------
      ! [-] Save statistics & results
      !-------------------------------------------------------------------------

      mphs%iinfo (IINFO_STRAT_FACTO )    = Aii_factors%choice
      mphs%iinfo (IINFO_STRAT_SCHUR )    = schur_strategy

      !  add schur
      Select Case( schur_strategy )
      Case(SCHUR_STRATEGY_Exact,&
           SCHUR_STRATEGY_ExactFromFacto)

         mphs%iinfo (IINFO_SCHUR_SIZEOF  )  = &
              Byte2MByte(XMPH_dm_sizeof( mphs%dm_schur ))
         mphs%iinfo (IINFO_SCHUR_NBENTRIES) = mphs%dm_schur%m  * mphs%dm_schur%n

      Case(SCHUR_STRATEGY_ApprxWithIluT  )

         mphs%iinfo (IINFO_SCHUR_SIZEOF  )  = &
              Byte2MByte(XMPH_sm_sizeof( mphs%sm_schur ))
         mphs%iinfo (IINFO_SCHUR_NBENTRIES) = mphs%sm_schur%nnz
     
      End Select

      ! add peak / memusage from sds 
      Call XMPH_sds_get_numstats (Aii_factors,piv)
      Call XMPH_sds_get_memstats (Aii_factors,sds_mem)
      Call MPH_mem_add2mem(mphs%mem,sds_mem)
      Call XMPH_sds_get_perfstats(Aii_factors, &
           mphs%rinfo (RINFO_FLOPS_FactoEstiElim),&
           mphs%rinfo (RINFO_FLOPS_FactoAssemb  ),&
           mphs%rinfo (RINFO_FLOPS_FactoElim    ))

      mphs%iinfo (IINFO_LCFACTO_NBPIVOTS) = INT(piv,4)
      mphs%iinfo (IINFO_LCFACTO_SIZEOF ) = &
           Byte2MByte(MPH_mem_getallusage(sds_mem))
      mphs%iinfo (IINFO_FAC_SDSMEMPEAK  ) = &
           Byte2MByte(MPH_mem_getallpeak   (sds_mem))

      !-------------------------------------------------------------------------
      ! [-] Finish
      !-------------------------------------------------------------------------

9999  Continue
      ! report status
      mphs%iinfo( IINFO_STATUS ) = iinfo
      Call XMPH_state_stepend(mphs)

    End Subroutine XMPH_maphys_Factorize

    ! [+] routine : MAPHYS_Precond ---------------------------------------------
    !
    !> Driver for the preconditioning step
    !!
    !! The preconditionner is a factorization of the assembled schur (or
    !! its approximation)
    !!
    !! @param[in,out ] mphs     the maphys instance 
    !!
    !! @author Yohan Lee-tin-yien
    !!
    Subroutine XMPH_maphys_Precond(mphs)
      !* Module(s) *!
      Use mph_log_mod
      Use XMPH_state_print_mod
      Use XMPH_state_mod
      Use XMPH_pcd_mod

      Implicit None

      !* Arguments *!
      Type(XMPH_maphys_t), Intent(inout) :: mphs

      !* Local variables *!

      ! Scalars
      Integer :: iinfo
      !- End of header ---------------------------------------------------------

      !-------------------------------------------------------------------------
      ! [-] Initialize variables
      !-------------------------------------------------------------------------

      ! Check previous status
      iinfo = MPH_SUCCESS
      CHCKASSRT( mphs%iinfo( IINFO_STATUS ) >= 0, iinfo)
      MPH_ONFAILURE_RETURN(iinfo)

      mphs%ikeep(IKEEP_CURJOB) = CURJOB_IsPrecond
      Call MPH_Time_start( mphs%rinfo (RINFO_TIMING_Precond) )

      Call XMPH_state_stepstart(mphs, mphs%ikeep, mphs%rkeep,iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)

      mphs%iinfo (IINFO_STRAT_PCD ) = mphs%ikeep(IKEEP_PCD_Strategy)

      !-------------------------------------------------------------------------
      ! [-] Form the preconditionner according to the selected strategy
      !-------------------------------------------------------------------------
      
      Select Case (mphs%ikeep(IKEEP_PCD_Strategy))
      Case (PCD_STRATEGY_isLocalExact )

         Call XMPH_pcd_LocalExact(mphs,iinfo)
         MPH_ONFAILURE_GOTO9999(iinfo)

      Case (PCD_STRATEGY_isLocalApprox )

         Call XMPH_pcd_LocalApprox(mphs,iinfo)
         MPH_ONFAILURE_GOTO9999(iinfo)

      Case (PCD_STRATEGY_isForcedByILUT ) 

         Call XMPH_pcd_LocalApproxFromApprox(mphs,iinfo)
         MPH_ONFAILURE_GOTO9999(iinfo)

      Case (PCD_STRATEGY_isNone)

         Call mph_log(MSG_Warning, "preconditioner deactivated")
         mphs%iinfo (IINFO_PCD_SDSMEMPEAK   ) = -1
         mphs%iinfo (IINFO_STRAT_PCDSDS  ) = -1
         mphs%iinfo (IINFO_PCD_ORDER     ) = -1
         mphs%iinfo (IINFO_PCD_SIZEOF   ) = -1
         mphs%iinfo (IINFO_PCD_NBPIVOTS  ) = -1
         mphs%iinfo (IINFO_PCD_PERKEPT   ) = -1

      Case Default

         Call mph_log(MSG_ERROR, "preconditioner unsupported")
         mphs%iinfo (IINFO_PCD_SDSMEMPEAK   ) = -1
         mphs%iinfo (IINFO_STRAT_PCDSDS  ) = -1
         mphs%iinfo (IINFO_PCD_ORDER     ) = -1
         mphs%iinfo (IINFO_PCD_SIZEOF   ) = -1
         mphs%iinfo (IINFO_PCD_NBPIVOTS  ) = -1
         mphs%iinfo (IINFO_PCD_PERKEPT   ) = -1
         
         CHCKASSRT( .False. , iinfo)
         MPH_ONFAILURE_GOTO9999(iinfo)

      End Select
      
      Call MPH_Time_stop( mphs%rinfo (RINFO_TIMING_Precond) )

      !-------------------------------------------------------------------------
      ! [-] Exit routine
      !-------------------------------------------------------------------------

9999  Continue

      mphs%iinfo( IINFO_STATUS ) = iinfo
      Call XMPH_state_stepend(mphs)

    End Subroutine XMPH_maphys_Precond


    ! [+] routine : MAPHYS_Solve -----------------------------------------------
    !
    !> Performs the solve step
    !!
    !! @param[in,out ] mphs          the maphys instance
    !!
    !! @author Azzam Haidar
    !! @author Yohan Lee-tin-yien
    !!
    Subroutine XMPH_maphys_Solve (mphs)

      !* Module(s) used *!

      Use mph_log_mod
      Use XMPH_schur_mod  
      Use XMPH_domainsls_mod
      Use XMPH_dense_matrix_mod
      Use XMPH_sparse_matrix_mod
      Use XMPH_state_mod
      Use XMPH_part_mod,        Only :  &
           XMPH_PART_DistRhs,    & ! routine(s)
           XMPH_PART_CollectSol
      Implicit None

      !* Arguments *!
      Type(XMPH_maphys_t), Intent(inout) :: mphs

      !* Local variables *!

      ! Scalars
      Integer                        :: iinfo  
      Integer                        :: iterativeSolver
      Integer                        :: master, rank
      Logical                        :: HasInitGuess
      Logical                        :: HasPcdFromPILUT
      Logical                        :: RHSIsCentralized

      ! Derived types
      Type(XMPH_dense_matrix_t)           :: bound_rhs
      Type(XMPH_dense_matrix_t)           :: bound_sol

      Type(XMPH_dense_matrix_t)           :: domain_rhs
      Type(XMPH_dense_matrix_t)           :: domain_sol

      Type(XMPH_dense_matrix_t)           :: gb_sol
      Type(XMPH_dense_matrix_t)           :: gb_rhs

      Type(XMPH_sparse_matrix_t)          :: gb_A

      !- End of header----------------------------------------------------------

      !-------------------------------------------------------------------------
      ! [-] Init
      !-------------------------------------------------------------------------

      ! Check previous status
      iinfo = MPH_SUCCESS
      CHCKASSRT( mphs%iinfo( IINFO_STATUS ) >= 0, iinfo)
      MPH_ONFAILURE_RETURN(iinfo)

      Call MPH_Time_reset(mphs%rinfo( RINFO_TIMING_SolveGathSol ))
      Call MPH_Time_reset(mphs%rinfo( RINFO_TIMING_SolveDistRHS ))
      Call MPH_Time_start(mphs%rinfo(RINFO_TIMING_Solve))
      
      mphs%ikeep(IKEEP_CURJOB) = CURJOB_IsSolve
      rank  =mphs%ikeep(IKEEP_MPIRANK)
      master=mphs%ikeep(IKEEP_HOSTRANK)


      ! read options
      Call XMPH_state_stepstart(mphs, mphs%ikeep, mphs%rkeep,iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)

      HasInitGuess = (mphs%IKEEP(IKEEP_ITS_InitGuess) == 1)
      RHSIsCentralized = (mphs%ikeep(IKEEP_INSYSTEM) == INSYSTEM_IsCentralized)
            
      ! set the iterative strategy
      IterativeSolver = mphs%ikeep(IKEEP_ITS_Solver)
      HasPcdFromPILUT = ( mphs%ikeep(IKEEP_SCHUR_STRATEGY) == &
           SCHUR_STRATEGY_ApprxWithIluT )

#if 0      
      If (HasPcdFromPILUT.And.(IterativeSolver == ITS_Solver_isPackCg))Then
         Call mph_log(MSG_WARNING, "Using GMRES instead of CG: &
              & CG does not support preconditioner from PILUT")
         IterativeSolver = ITS_Solver_isPackGMRES
         mphs%ikeep(IKEEP_ITS_Solver) = IterativeSolver
      End If
#endif
      ! check a few array association
      If (RHSIsCentralized.And.(rank == master))Then
         CHCKASSRT( Associated(mphs%rhs) , iinfo )
         If (iinfo < 0 ) Goto 9999
         CHCKASSRT( Associated(mphs%sol) , iinfo )
         If (iinfo < 0 ) Goto 9999
      Endif      

      ! nullify local types
      Call XMPH_dm_Nullify( bound_rhs , iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)

      Call XMPH_dm_Nullify( bound_sol , iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)
      
      Call XMPH_dm_Nullify( domain_rhs, iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)
      
      Call XMPH_dm_Nullify( domain_sol, iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)

      Call XMPH_dm_Nullify( gb_sol    , iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)

      Call XMPH_dm_Nullify( gb_rhs    , iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)


      !-------------------------------------------------------------------------
      ! [-] Generate the right-hand-sides & the init guesses
      !-------------------------------------------------------------------------

      !---------------------------------------------------------------
      ! [--] Generate the domain rhs 
      !---------------------------------------------------------------

      If (RHSIsCentralized)Then

         Call MPH_Time_start( mphs%rinfo( RINFO_TIMING_SolveDistRHS ))
         Call XMPH_PART_distRhs(             & ! intents 
              mphs%part_rhs, mphs%comm, & ! in
              mphs%lc_domain, mphs%rhs, & !    
              domain_rhs, iinfo         & ! out
              )
         Call MPH_Time_stop( mphs%rinfo( RINFO_TIMING_SolveDistRHS ))

         MPH_ONFAILURE_GOTO9999(iinfo)
      Else

         Call XMPH_DM_Create &
              (domain_rhs, mphs%lc_domain%myndof, &
               1,mphs%lc_domain%myndof,iinfo)
         MPH_ONFAILURE_GOTO9999(iinfo)
         Call XMPH_xcopy(domain_rhs%v,mphs%rhs, domain_rhs%m)

      EndIf

      If (RHSIsCentralized .And. HasInitGuess)Then
         Call XMPH_PART_distRhs(             & ! intents 
           mphs%part_rhs, mphs%comm, & ! in
           mphs%lc_domain, mphs%sol, & !     
           domain_sol, iinfo         & ! out
           )
         MPH_ONFAILURE_GOTO9999(iinfo)
      End If


      !---------------------------------------------------------------
      ! [--] Generate the rhs for the interface linear system 
      !---------------------------------------------------------------
      !
      ! bound_rhs := SUM [ domain_RHSb - (Abi . Aii^-1 . domain_RHSi ) ]
      ! or 
      ! bound_rhs := domain_RHSb - SUM [ (Abi . Aii^-1 . domain_RHSi ) ]
      ! 
      ! depending on how the global RHS was distributed.
      !
      Call MPH_Time_start( mphs%rinfo( RINFO_TIMING_SolveGenRHS))
      Call XMPH_schur_generate_RHS( & ! intents
           mphs,                & ! inout
           domain_rhs,          & ! in
           bound_rhs,           & ! out
           iinfo                &
           )
      Call MPH_Time_stop( mphs%rinfo( RINFO_TIMING_SolveGenRHS))
      MPH_ONFAILURE_GOTO9999(iinfo)

      !---------------------------------------------------------------
      ! [--] Extract the initial guess on the interface linear system 
      !---------------------------------------------------------------

      If (HasInitGuess)Then
         
         Call XMPH_schur_Generate_InitGuess(mphs,domain_sol,bound_sol,iinfo)
         MPH_ONFAILURE_GOTO9999(iinfo) 

      End If

      !-------------------------------------------------------------------------
      ! [-] Solve the system on the interfaces (iterative method)
      !-------------------------------------------------------------------------

      Call XMPH_SCHUR_InitTimers(mphs)
      Call MPH_Time_start( mphs%rinfo( RINFO_TIMING_SolveITS) )

      Select Case ( IterativeSolver )
      Case (ITS_Solver_isPackCg   )
         Call XMPH_schur_CG_Solve (mphs, bound_rhs,bound_sol, iinfo)
      Case (ITS_SOLVER_iSPackGmres)
         Call XMPH_schur_GMRES_Solve(mphs,bound_rhs,bound_sol, iinfo)
      Case (ITS_SOLVER_iSPackFGmres)
         Call XMPH_schur_FGMRES_Solve(mphs,bound_rhs,bound_sol, iinfo)
      Case Default 
         iinfo = -3
      End Select
      
      Call MPH_Time_stop( mphs%rinfo( RINFO_TIMING_SolveITS) )

      MPH_ONFAILURE_GOTO9999(iinfo)

      !-------------------------------------------------------------------------
      ! [-] Solve the system on the interior  (direct method)
      ! [-] knowing the solution the boundary
      !-------------------------------------------------------------------------

      Call MPH_Time_start( mphs%rinfo( RINFO_TIMING_SolveSDS) )
      Call XMPH_domainsls_Solve_interior( & !
           mphs%sls_domain%sds,   & !
           mphs%sm_Aib,           & !
           domain_rhs, bound_sol, & !
           domain_sol, iinfo)
      Call MPH_Time_stop( mphs%rinfo( RINFO_TIMING_SolveSDS) )

      MPH_ONFAILURE_GOTO9999(iinfo)

      !-------------------------------------------------------------------------
      ! [-] Gather the solution & compute related statistics
      !-------------------------------------------------------------------------
      
      !---------------------------------------------------------------
      ! [--] Gather the solution 
      !---------------------------------------------------------------

      If (RHSIsCentralized)Then

         Call MPH_Time_start(mphs%rinfo( RINFO_TIMING_SolveGathSol ))
         Call XMPH_PART_collectSol(       & ! intents
              mphs%part_rhs,              & ! inout
              mphs%comm,  mphs%lc_domain, & ! in
              domain_sol,                 &              
              gb_sol, iinfo               & ! out
              )
         MPH_ONFAILURE_GOTO9999(iinfo) 
         Call MPH_Time_stop(mphs%rinfo( RINFO_TIMING_SolveGathSol ))

      End If

      Call MPH_Time_stop(mphs%rinfo(RINFO_TIMING_Solve))

      !---------------------------------------------------------------
      ! [--] Compute statistics
      !---------------------------------------------------------------

      !
      mphs%iinfo (IINFO_STRAT_ITS ) = IterativeSolver

      mphs%rinfog(RINFOG_BckwrdErrorDist) = XMPH_EstimateError &
           ( mphs%sls_domain%sm_A, domain_sol,domain_rhs, mphs )

      If (RHSIsCentralized.And.(rank == master)) Then

         ! gb_A
         Call XMPH_sm_ijv                        & ! intents
              (gb_A,                             & ! out
              mphs%rows, mphs%cols, mphs%values, & ! in
              mphs%n, mphs%n,mphs%nnz, mphs%sym, & ! 
              iinfo                              & ! out           
              )
         MPH_ONFAILURE_GOTO9999(iinfo)

         ! gb_rhs
         gb_rhs%m  = gb_A%m
         gb_rhs%n  = 1
         gb_rhs%ld = gb_A%m
         gb_rhs%v  => mphs%rhs(1:gb_rhs%m)

         mphs%rinfog(RINFOG_BACKWARD_ERROR) = &
              XMPH_EstimateErrorCent(gb_A,gb_sol,gb_rhs)

      End If

      ! Save total timings

      Call MPH_Time_stop(mphs%rinfo(RINFO_TIMING_Total))

      mphs%rinfo (RINFO_TIMING_AllSteps) = &
           mphs%rinfo (RINFO_TIMING_Analysis) + &
           mphs%rinfo (RINFO_TIMING_Facto   ) + &
           mphs%rinfo (RINFO_TIMING_Precond ) + &
           mphs%rinfo (RINFO_TIMING_Solve   )

      !-------------------------------------------------------------------------
      ! [-] Finish
      !-------------------------------------------------------------------------

      ! Save solution
      If (RHSIsCentralized .And.(rank == master))&
         Call XMPH_xcopy(mphs%sol,gb_sol%v, gb_sol%m)

9999  Continue

      ! Free memory
      Call XMPH_dm_Free(bound_rhs, iinfo)
      Call XMPH_dm_Free(bound_sol, iinfo)
      Call XMPH_dm_Free(domain_rhs, iinfo)
      Call XMPH_dm_Free(domain_sol, iinfo)
      If (rank == master ) Call XMPH_dm_Free(gb_sol, iinfo)

      ! report status
      mphs%iinfo( IINFO_STATUS ) = iinfo
      Call XMPH_state_stepend(mphs)

    end Subroutine XMPH_maphys_Solve


    ! [+] routine : MAPHYS_Exit  -----------------------------------------------
    !
    !> Destroy a maphys_instance
    !! 
    !! @param[in,out ] mphs     the maphys instance to finalize
    !!
    !! @author Yohan Lee-tin-yien
    !!
    Subroutine XMPH_maphys_Exit(mphs) ! inout

      !* Modules *!
      Use XMPH_maphys_type
      Use XMPH_maphys_aux_mod
      Implicit None

      !* Arguments *!
      Type(XMPH_maphys_t), intent(inout) :: mphs

      !* Local variables *!

      ! Scalars
      Integer :: iinfo 

      !- End of header ---------------------------------------------------------

      !-------------------------------------------------------------------------
      ! [-] Free everything
      !-------------------------------------------------------------------------

      Call XMPH_maphys_free(mphs, iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)

      !-------------------------------------------------------------------------
      ! [-] Exit routine
      !-------------------------------------------------------------------------

9999  Continue
      ! report status
      mphs%iinfo( IINFO_STATUS ) = iinfo

    End Subroutine XMPH_maphys_Exit

    !> Estimate the error on the solution
    !!
    !! Compute | A.x - b |/|b| (with |.| norm 2)
    !! where A is the initial matrix, which was distributed.
    !! @Return the norm2 of the error, on failure return "-1.d0"
    !! 
    !! @param[in] sm  the local matrix 
    !! @param[in] sol the local part of the solution
    !! @param[in] rhs the local part of right-hand-side 
    !! @param[in,out] mphs the maphys instance, uses components :
    !!       - [in] comm the MPI communictor on the domains
    !!       - [in] lc_domain the domain description
    !!       - [in,out] intrfbuff the buffer for communications on the interface
    !! 
    Real*8 Function XMPH_EstimateError &
         ( sm, sol, rhs, mphs)
      
      !* Module(s) *!

      Use XMPH_sparse_matrix_mod
      Use XMPH_dense_matrix_mod
      Use XMPH_part_mod
      Use XMPH_schur_aux_mod
      Implicit None

      !* Argument(s) *!

      Type(XMPH_sparse_matrix_t), Intent(in) :: sm
      Type(XMPH_dense_matrix_t), Intent(in) :: sol
      Type(XMPH_dense_matrix_t), Intent(in) :: rhs
      Type(XMPH_maphys_t), Intent(inout) :: mphs

      !* External routines *!
      Include 'mpif.h'
      XMPH_FLOAT, External :: XMPH_DOT
      
      !* Local variable(s) *!
      Integer :: iinfo, comm
      MPH_INT :: i
      MPH_INT :: lim
      Real*8  :: squareOfNormRHS
      Real*8  :: squareOfUpper
      Real*8  :: squareOfNormRHSRecv
      Real*8  :: squareOfUpperRecv
      Type(XMPH_dense_matrix_t) :: true_rhs
      Type(XMPH_dense_matrix_t) :: dm_dom ! temporary vector on the domain
      Type(XMPH_dense_matrix_t) :: dm_itf ! temporary vector on the interface


      !- End of header----------------------------------------------------------

      !-------------------------------------------------------------------------
      ! [-] Init
      !-------------------------------------------------------------------------

      ! Init result to the error value
      XMPH_EstimateError = -1.d0

      ! Limitations of current routine
      ASSRT(rhs%n==1)

      ! Init 
      iinfo = MPH_SUCCESS
      lim = mphs%lc_domain%myndofinterior+1
      comm = mphs%comm
      
      Call XMPH_DM_Nullify(true_rhs,iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)

      Call XMPH_DM_Nullify(dm_dom  ,iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)

      Call XMPH_DM_Nullify(dm_itf  ,iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)

      squareOfUpper = 0
      squareOfnormRHSRecv = 0
      !-------------------------------------------------------------------------
      ! [-] compute |b|^2 
      !-------------------------------------------------------------------------
      !> @note 
      !! This step is validated, while weight(i) == 0 or 1.

      Call XMPH_dm_dup(rhs,true_rhs,iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)

      dm_itf%m  = mphs%lc_domain%myndofintrf
      dm_itf%n  = true_rhs%n
      dm_itf%ld = true_rhs%ld
      dm_itf%v => true_rhs%v &
           (mphs%lc_domain%myndofinterior+1 : true_rhs%m*true_rhs%ld)
      Do i=1,dm_itf%m
         dm_itf%v(i) = mphs%lc_domain%weight(i) * dm_itf%v(i)
      End Do
      squareOfNormRHS = Real(XMPH_DOT(&
           true_rhs%m, true_rhs%v(1), 1, true_rhs%v(1),1 &
           ),8)
      Call MPI_AllReduce(squareOfNormRHS,squareOfNormRHSRecv, &
           1, MPI_DOUBLE_PRECISION,MPI_SUM,comm,iinfo)
      ASSRT(iinfo == MPI_SUCCESS)
      squareOfNormRHS = squareOfNormRHSRecv
      !-------------------------------------------------------------------------
      ! [-] dm_dom <-- A.x - b
      !-------------------------------------------------------------------------

      Call XMPH_DM_Create(dm_dom,sol%m,1,sol%ld,iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)

      Call XMPH_sm_VectorProduct &
           (sm, 0, 0, sol, dm_dom, iinfo )
      MPH_ONFAILURE_GOTO9999(iinfo)
     
      Do i=1,dm_dom%m
         dm_dom%v(i) = dm_dom%v(i) - rhs%v(i)
      End Do

      !-------------------------------------------------------------------------
      ! [-] Handle the interface
      !-------------------------------------------------------------------------

      Call XMPH_DM_Nullify(dm_itf,iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)

      dm_itf%m  = mphs%lc_domain%myndofintrf
      dm_itf%n  = dm_dom%n
      dm_itf%ld = dm_dom%ld
      dm_itf%v => dm_dom%v(mphs%lc_domain%myndofinterior+1:dm_dom%m)

      Call XMPH_SCHUR_UpdateVector(mphs,dm_itf,iinfo)
      MPH_ONFAILURE_GOTO9999(iinfo)

      Do i=1,dm_itf%m
         dm_itf%v(i) = dm_itf%v(i) * mphs%lc_domain%weight(i)
      End Do

      !-------------------------------------------------------------------------
      ! [-] squareOfUpper <-- |A.x - b|^2
      !-------------------------------------------------------------------------

      squareOfUpper = Real(XMPH_DOT(dm_dom%m,dm_dom%v(1), 1, dm_dom%v(1),1 ),8)
      Call MPI_AllReduce(squareOfUpper,squareOfUpperRecv,  & 
           1, MPI_DOUBLE_PRECISION,MPI_SUM,comm,iinfo)
      ASSRT(iinfo == MPI_SUCCESS)
      squareOfUpper = squareOfUpperRecv
      !-------------------------------------------------------------------------
      ! [-] Get the error
      !-------------------------------------------------------------------------

      XMPH_EstimateError = Sqrt(squareOfUpper) / Sqrt(squareOfNormRHS)

      !-------------------------------------------------------------------------
      ! [-] Exit
      !-------------------------------------------------------------------------
      
9999  Continue
      Call XMPH_DM_Nullify(dm_itf,iinfo)
      Call XMPH_DM_Free(dm_dom,iinfo)
      Call XMPH_DM_Free(true_rhs,iinfo)
      
    End Function XMPH_EstimateError

    !> Estimate the error on the solution with centralized data
    !!
    !! Compute | A.x - b |/|b| (with |.| norm 2)
    !! where A is the initial matrix, which was distributed.
    !! @Return the norm2 of the error, on failure return "-1.d0"
    !! 
    !! @param[in] sm  the global matrix 
    !! @param[in] sol the global solution
    !! @param[in] rhs the global right-hand-side 
    !!
    Real*8 Function XMPH_EstimateErrorCent &
         ( sm, sol, rhs )

      !* Module(s) *!

      Use XMPH_sparse_matrix_mod
      Use XMPH_dense_matrix_mod

      Implicit None

      !* Argument(s)* !

      Type(XMPH_sparse_matrix_t), Intent(in) :: sm
      Type(XMPH_dense_matrix_t), Intent(in)  :: sol
      Type(XMPH_dense_matrix_t), Intent(in)  :: rhs

      !* External functions *!
      XMPH_FLOAT, External :: XMPH_DOT

      !* Local variables *!
      MPH_INT :: i,n
      Integer :: iinfo
      Type(XMPH_dense_matrix_t) :: tmp
      Real*8  :: squareOfNormRHS
      Real*8  :: squareOfUpper

      !- End of header----------------------------------------------------------

      ! Init

      XMPH_EstimateErrorCent = -1.d0
      n = sol%m
      Call XMPH_dm_Nullify(tmp,iinfo )
      MPH_ONFAILURE_GOTO9999(iinfo)

      Call XMPH_dm_Create(tmp,n,1,n,iinfo )
      MPH_ONFAILURE_GOTO9999(iinfo)

      ! tmp <-- A.x

      Call XMPH_sm_VectorProduct(sm, 0, 0, sol, tmp, iinfo )
      MPH_ONFAILURE_GOTO9999(iinfo)
      
      ! tmp <-- A.x - b
      Do i=1, n
         tmp%v(i) = tmp%v(i) - rhs%v(i)
      End Do
      
      ! | A.x - b |/|b|

      squareOfUpper   = Real(XMPH_DOT(n, tmp%v(1), 1, tmp%v(1),1 ),8)
      squareOfNormRHS = Real(XMPH_DOT(n, rhs%v(1), 1, rhs%v(1),1 ),8)
      XMPH_EstimateErrorCent = SQRT( squareOfUpper ) / SQRT( squareOfNormRHS ) 
      
      ! Exit
9999  Continue
      Call XMPH_dm_free(tmp,iinfo)

    End Function XMPH_EstimateErrorCent

    
    
  End Module XMPH_maphys_mod

