! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"

! [+] module : XMPH_maphys_aux_mod ---------------------------------------------
!
!> Auxialiary module for maphys
!!
!! Basically contains useful routines for maphys.
!!
Module XMPH_maphys_aux_mod
  Implicit None

  !* Private constants *!
  Character(len=MAPHYS_STRL), Private, Parameter :: FLNAME= &
       "XMPH_ARITHmph_maphys_aux_mod.F90"

  ! List of routines
  Public :: XMPH_MAPHYS_Nullify
  Public :: XMPH_MAPHYS_Free

  Public :: XMPH_MAPHYS_Set_default_icntl
  Public :: XMPH_MAPHYS_Set_default_rcntl

  Public :: XMPH_MAPHYS_PretreatInputMatrix
  
Contains

  ! [+] routine : XMPH_MAPHYS_Nullify  ----------------------------------------------
  !
  !> nullify a maphys_instance
  !! 
  !! @param[in,out ] mphs     the maphys instance to nullify
  !!
  !! @author Yohan Lee-tin-yien
  !! @todo   nullify derived type : mphs%dls_precond_schur%dds
  !!
  Subroutine XMPH_maphys_Nullify & ! intents
       (mphs ,           & ! inout
       info )              ! out
    
    !* Modules *!

    Use XMPH_maphys_type
    Use XMPH_sparse_matrix_mod, Only : &
         XMPH_sm_nullify ! routine
    Use XMPH_dense_matrix_mod, Only : &
         XMPH_dm_nullify ! routine

  
    implicit none

    !* Subroutine arguments *!
    Type(XMPH_maphys_t), intent(inout) :: mphs
    integer       , intent(  out) :: info

    !* Local variables *!

    ! contants
    integer       , parameter :: DEFAULT_INT   = -9999
    real(kind=8)  , parameter :: DEFAULT_REAL  = -9999.d0

    ! others
    integer                       :: iinfo

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] USER PARAMETERS
    !---------------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! [1.1] MPI communicator (do not touch)
    !-------------------------------------------------------------------
    ! mphs%comm = 

    !-------------------------------------------------------------------
    ! [1.2] Input matrix (in coordinate format)
    !-------------------------------------------------------------------
    mphs%sym = DEFAULT_INT
    mphs%n   = DEFAULT_INT
    mphs%nnz = DEFAULT_INT
    nullify( mphs%rows )
    nullify( mphs%cols )
    nullify( mphs%values )
    mphs%write_matrix = ""

    !-------------------------------------------------------------------
    ! [1.3] Right-hand-side & Solution (in dense format ordered by columns)
    !-------------------------------------------------------------------

    mphs%nrhs  = DEFAULT_INT
    nullify( mphs%rhs ) 
    nullify( mphs%sol ) 

    !-------------------------------------------------------------------
    ! [1.4] Controls
    !-------------------------------------------------------------------

    mphs%job    = DEFAULT_INT
    mphs%icntl  = DEFAULT_INT
    mphs%rcntl  = DEFAULT_REAL

    !-------------------------------------------------------------------
    ! [1.5] Statistics
    !-------------------------------------------------------------------

    ! on this process (MPI)
    mphs%iinfo = -1
    mphs%rinfo = -1

    ! on all process (MPI)
    mphs%iinfog = -1
    mphs%rinfog = -1

    mphs%iinfomin = -1
    mphs%iinfomax = -1
    mphs%iinfoavg = -1
    mphs%iinfosig = -1
    mphs%rinfomin = -1
    mphs%rinfomax = -1
    mphs%rinfoavg = -1
    mphs%rinfosig = -1

    !---------------------------------------------------------------------------
    ! [2] Internal working data
    !---------------------------------------------------------------------------


    !-------------------------------------------------------------------
    ! [2.1] internal controls
    !-------------------------------------------------------------------

    mphs%ikeep  = DEFAULT_INT
    mphs%rkeep  = DEFAULT_REAL

    !-------------------------------------------------------------------
    ! [2.2]  Description of the Non-overlapping domain Decomposition 
    !-------------------------------------------------------------------

    !---------------------------------------------------------
    ! [2.2.1] Local domain description (interface + interior)
    !---------------------------------------------------------
    ! Scalars
    mphs%lc_domain%lenindintrf      = DEFAULT_INT
    mphs%lc_domain%myint_Lg         = DEFAULT_INT
    mphs%lc_domain%mynbvi           = DEFAULT_INT
    mphs%lc_domain%myndof           = DEFAULT_INT
    mphs%lc_domain%myndofinterior   = DEFAULT_INT
    mphs%lc_domain%myndofintrf      = DEFAULT_INT
    mphs%lc_domain%myndoflogicintrf = DEFAULT_INT
    mphs%lc_domain%mysizeIntrf      = DEFAULT_INT

    ! Arrays

    Nullify( mphs%lc_domain%myindexVi    )
    Nullify( mphs%lc_domain%myindexintrf )
    Nullify( mphs%lc_domain%myinterface  )
    Nullify( mphs%lc_domain%mylogicintrf )
    Nullify( mphs%lc_domain%myptrindexVi )
    Nullify( mphs%lc_domain%weight       )

    !---------------------------------------------------------
    ! [2.2.2] Necessary data to part the right hand side
    !---------------------------------------------------------
    ! Scalars

    mphs%part_rhs%gballintrf = DEFAULT_INT
    mphs%part_rhs%gballndof = DEFAULT_INT

    mphs%part_rhs%combivtxpcsz = DEFAULT_INT
    mphs%part_rhs%rhsway       = DEFAULT_INT

    ! Arrays
    Nullify( mphs%part_rhs%domLg            )
    Nullify( mphs%part_rhs%domintdof        )
    Nullify( mphs%part_rhs%domintrfdof      )
    Nullify( mphs%part_rhs%domlogicintrfdof )
    Nullify( mphs%part_rhs%metperm          )
    Nullify( mphs%part_rhs%procintdisp      )
    Nullify( mphs%part_rhs%procintrfdisp    )
    Nullify( mphs%part_rhs%scatindices      )
    Nullify( mphs%part_rhs%scatlogicindices )
    Nullify( mphs%part_rhs%gbtoloc      )

    !---------------------------------------------------------
    ! [2.2.3] Blocs on the matrix of the domain 
    !---------------------------------------------------------

    ! [  Aii     Aib ] 
    ! [  Abi     Abb ] 
    ! - i: interior  related (nodes inside the domain)
    ! - b: border    related (nodes on its interface )

    Call XMPH_sm_nullify(mphs%sm_Aii,iinfo)
    Call XMPH_sm_nullify(mphs%sm_Aib,iinfo)
    Call XMPH_sm_nullify(mphs%sm_Abi,iinfo)
    Call XMPH_sm_nullify(mphs%sm_Abb,iinfo)

    !---------------------------------------------------------
    ! [2.2.4] saved permutation / inverse permutation
    !---------------------------------------------------------
    !lty> deprecated -- contained in part_rhs
    ! nullify( mphs%domain_perm )
    ! nullify( mphs%domain_invperm )
    !lty<

    !---------------------------------------------------------
    ! [2.2.5] saved scaling
    !---------------------------------------------------------
    nullify( mphs%row_scaling )  
    nullify( mphs%col_scaling )  

    !-------------------------------------------------------------------
    ! [2.3] Linear system on a domain
    !-------------------------------------------------------------------

    !---------------------------------------------------------
    ! [2.3.1] Its interior is solved with a sparse direct solver
    !---------------------------------------------------------

    ! Scalars 
    mphs%sls_domain%rhs_is_sparse = MPH_FALSE

    ! Arrays
    Nullify(mphs%sls_domain%sm_A  )
    Nullify(mphs%sls_domain%sds   )
    Nullify(mphs%sls_domain%dm_rhs)
    Nullify(mphs%sls_domain%sm_rhs)
    Nullify(mphs%sls_domain%dm_sol)

    !-------------------------------------------------------------------
    ! [2.3.2] On the interface : the schur complement linear system
    ! (solved with an iterative method)
    !-------------------------------------------------------------------

    ! local schur complements matrices (exact or approximate)
    Call XMPH_dm_nullify(mphs%dm_schur,iinfo)
    Call XMPH_sm_nullify(mphs%sm_schur,iinfo)

    !** preconditioners for the iterative method **!
    !* dense preconditioner *!

    ! Arrays
    Nullify(mphs%dls_precond_schur%IPIV)

    ! Derived Types
    !@ Type(XMPH_dense_direct_solver_t)   :: dds 
    Call XMPH_dm_nullify(mphs%dls_precond_schur%dm_A,iinfo)
    Call XMPH_dm_nullify(mphs%dls_precond_schur%dm_B,iinfo)
    Call XMPH_dm_nullify(mphs%dls_precond_schur%dm_X,iinfo)

    !* sparse preconditioner *!

    ! Scalars 
    mphs%sls_precond_schur%rhs_is_sparse = MPH_FALSE

    ! Arrays
    Nullify(mphs%sls_precond_schur%sm_A  )
    Nullify(mphs%sls_precond_schur%sds   )
    Nullify(mphs%sls_precond_schur%dm_rhs)
    Nullify(mphs%sls_precond_schur%sm_rhs)
    Nullify(mphs%sls_precond_schur%dm_sol)

    !-------------------------------------------------------------------
    ! [3] Finish
    !-------------------------------------------------------------------
    
    ! force the success of the routine
    info          = 0
    mphs%iinfo(1) = 0


  end Subroutine XMPH_maphys_Nullify

! [+] routine : XMPH_maphys_Free  ---------------------------------------------------
!
!> Free a maphys_instance
!! 
!! @param[in,out ] mphs     the maphys instance to nullify
!!
!! @author Yohan Lee-tin-yien
!! @todo   Implement
Subroutine XMPH_maphys_Free    & ! intents
     (mphs ,           & ! inout
     info )              ! out

  Use XMPH_maphys_type
  Use MPH_maphys_enum
  Use XMPH_sparse_matrix_mod, Only : &
       XMPH_sm_nullify,   &
       XMPH_sm_free
  Use XMPH_dense_matrix_mod, Only : &
       XMPH_dm_nullify,   &
       XMPH_dm_free
  Use XMPH_sds_mod, Only : &
       XMPH_sds_finalize
  Use XMPH_dds_mod, Only : &
       XMPH_dds_exit
  Use XMPH_sls_mod, Only : &
       XMPH_sls_free

  Implicit None
  Include 'mpif.h'

  !* Subroutine arguments *!
  Type(XMPH_maphys_t), Intent(inout) :: mphs
  Integer       , Intent(  out) :: info

  !* Local variables *!

  ! others
  Integer                       :: iinfo

  !- End of header -------------------------------------------------------------


  !-----------------------------------------------------------------------------
  ! [-] USER PARAMETERS (nothing to do)
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! [-] Internal working data
  !-----------------------------------------------------------------------------

  !-------------------------------------------------------------------
  ! [--]  Description of the Non-overlapping domain Decomposition 
  !-------------------------------------------------------------------

  !---------------------------------------------------------
  ! [---] Local domain description (interface + interior)
  !---------------------------------------------------------
  !@ Type(XMPH_maphys_domain_t)     :: lc_domain 

  If (Associated(mphs%lc_domain%weight))&
       Deallocate( mphs%lc_domain%weight)

  If (Associated(mphs%lc_domain%myindexintrf))&
       Deallocate( mphs%lc_domain%myindexintrf)

  If (Associated(mphs%lc_domain%myinterface))&
       Deallocate( mphs%lc_domain%myinterface )

  If (Associated(mphs%lc_domain%mylogicintrf))&
       Deallocate( mphs%lc_domain%mylogicintrf )

  If (Associated(mphs%lc_domain%myindexVi))&
       Deallocate( mphs%lc_domain%myindexVi )

  If (Associated(mphs%lc_domain%myptrindexVi))&
       Deallocate( mphs%lc_domain%myptrindexVi )


  !---------------------------------------------------------
  ! [---] Necessary data to part the right hand side
  !---------------------------------------------------------
  !@ Type(XMPH_maphys_rhs_partition_t) :: part_rhs
  ! @warning we ignore errors here
  If (Associated (mphs%part_rhs%procintdisp )) & 
       Deallocate(mphs%part_rhs%procintdisp      )

  If (Associated(mphs%part_rhs%procintrfdisp)) &
       Deallocate(mphs%part_rhs%procintrfdisp    )

  If (Associated(mphs%part_rhs%domLg)) &
       Deallocate(mphs%part_rhs%domLg            )

  If (Associated(mphs%part_rhs%domintdof)) &
       Deallocate(mphs%part_rhs%domintdof        )

  If (Associated(mphs%part_rhs%domintrfdof)) &
      Deallocate(mphs%part_rhs%domintrfdof       )

  If (Associated(mphs%part_rhs%metperm)) &
       Deallocate(mphs%part_rhs%metperm          )

  If (Associated(mphs%part_rhs%scatindices)) &
       Deallocate(mphs%part_rhs%scatindices      )

  If (Associated(mphs%part_rhs%scatlogicindices)) &
       Deallocate(mphs%part_rhs%scatlogicindices )

  If (Associated(mphs%part_rhs%domlogicintrfdof)) &
       Deallocate(mphs%part_rhs%domlogicintrfdof )

  If (Associated(mphs%part_rhs%gbtoloc))&
       Deallocate( mphs%part_rhs%gbtoloc     )


  !---------------------------------------------------------
  ! [---] Blocs on the matrix of the domain 
  !---------------------------------------------------------

  ! [  Aii     Aib ] 
  ! [  Abi     Abb ] 
  ! - i: interior  related (nodes inside the domain)
  ! - b: bound     related (nodes on its interface )

  Call XMPH_sm_free(mphs%sm_Aii,iinfo)
  Call XMPH_sm_free(mphs%sm_Aib,iinfo)
  Call XMPH_sm_free(mphs%sm_Abi,iinfo)
  Call XMPH_sm_free(mphs%sm_Abb,iinfo)

  !---------------------------------------------------------
  ! [---] saved scaling
  !---------------------------------------------------------

  If (Associated( mphs%row_scaling )) Deallocate( mphs%row_scaling )  
  If (Associated( mphs%col_scaling )) Deallocate( mphs%col_scaling )  

  !-------------------------------------------------------------------
  ! [--] Linear system on a domain
  !-------------------------------------------------------------------

  !---------------------------------------------------------
  ! [---] Its interior is solved with a sparse direct solver
  !---------------------------------------------------------

  !@ Type(XMPH_sparse_linear_system_t) :: sls_domain  
  Call XMPH_SLS_Free( mphs%sls_domain , iinfo )

  !-------------------------------------------------------------------
  ! [---] On the interface : the schur complement linear system
  ! (solved with an iterative method)
  !-------------------------------------------------------------------

  ! local schur complements matrices (exact or approximate)
  ! (schur allocated by sds)
  Call XMPH_dm_nullify(mphs%dm_schur,iinfo)
  Call XMPH_sm_free(mphs%sm_schur,iinfo)

  ! preconditioners for the iterative method

  Select Case (mphs%ikeep(IKEEP_Pcd_Strategy))
  Case(PCD_STRATEGY_isLocalExact)
     ! dls_precond_schur

     If ( mphs%ikeep(IKEEP_ITS_MatVect) == MAT_VECT_isImplicit ) Then
        ! memory already freed by call to XMPH_SLS_Free(mphs%sds)
        Call XMPH_dm_nullify(mphs%dls_precond_schur%dm_A,iinfo)
     Else
        Call XMPH_dm_free(mphs%dls_precond_schur%dm_A,iinfo)
     End If

     Call XMPH_dm_free(mphs%dls_precond_schur%dm_B,iinfo)
     Call XMPH_dm_free(mphs%dls_precond_schur%dm_X,iinfo)
     Call XMPH_DDS_Exit(mphs%dls_precond_schur%dds,iinfo)
     If (Associated(mphs%dls_precond_schur%ipiv)) &
          Deallocate( mphs%dls_precond_schur%ipiv)
  Case(PCD_STRATEGY_isLocalApprox)
     Call XMPH_SLS_Free( mphs%sls_precond_schur , iinfo )
  End Select

  !------------------------------------------------------------------
  ! [-] Finish
  !------------------------------------------------------------------

  ! This routine always return SUCCESS ( Do not check for errors )
  info = 0

End Subroutine XMPH_maphys_free


  ! [+] routine : XMPH_maphys_Set_default_icntl  -------------------------------
  !
  !> initialize controls (integers)
  !!
  !! give default values to ICNTLs
  !!
  !! @param[in,out ] icntl   integer controls
  !! @param[   out ] info    routine status
  !!
  !! @author Yohan Lee-tin-yien
  !! @todo Complete 
  Subroutine XMPH_maphys_Set_default_icntl(icntl,info)

    Use MPH_maphys_enum ! only ICNTL_*
    Use XMPH_sds_mod, Only : &
        XMPH_sds_t, & ! structure 
        XMPH_sds_select
    Implicit None

    !* Subroutine arguments *!

    Integer, Intent(inout) :: icntl( MAPHYS_ICNTL_SIZE )
    Integer, Intent(  out) :: info

    !* Local Variables *!

    ! Constants
    Integer, Parameter :: ZERO = 0

    ! Scalars
    Integer :: iinfo
    Integer :: i

    ! Derived types
    Type(XMPH_sds_t) :: sds_dummy

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [0] Initialize data
    !---------------------------------------------------------------------------
    !

    ! Initialize by default to ZERO
    Do i = 1, MAPHYS_ICNTL_SIZE
       icntl(i) = ZERO
    End Do

    !---------------------------------------------------------------------------
    ! [1] Set static values
    !---------------------------------------------------------------------------
    !

    ! outputs
    icntl( ICNTL_OUTPUT_ErrUnit   ) = 0
    icntl( ICNTL_OUTPUT_WarnUnit  ) = 0
    icntl( ICNTL_OUTPUT_StdUnit   ) = 6
    icntl( ICNTL_OUTPUT_Verbosity ) = 3

    icntl( ICNTL_PRINT_Cntls      ) = PRINT_Not
    icntl( ICNTL_PRINT_Infos      ) = PRINT_Not

    ! shift
    icntl( ICNTL_SHIFT_Do ) = MPH_FALSE

    ! partionner strategy
    icntl( ICNTL_PART_Strategy ) = 2

    ! Iterative Solver
    icntl( ICNTL_ITS_Solver      ) = 3
    icntl( ICNTL_ITS_OrtStrategy ) = 3
    icntl( ICNTL_ITS_InitGuess   ) = 0
    icntl( ICNTL_ITS_MaxIters    ) = 0
    icntl( ICNTL_ITS_ResStrategy ) = 0 
    icntl( ICNTL_ITS_Restart     ) = 0

    icntl( ICNTL_ITS_MatVect     ) = 0

    ! Schur Complement 
    icntl( ICNTL_SCHUR_Strategy    ) = 0

    icntl( ICNTL_EXSCHUR_BlocSize  ) = 50

    icntl( ICNTL_ILU_LUFillPct    ) = 0
    icntl( ICNTL_ILU_SCHURFillPct ) = 0

    icntl( ICNTL_ILU_LUFill     ) = -1
    icntl( ICNTL_ILU_SCHURFill  ) = -1

    ! 2nd lvl of parallelism
    icntl( ICNTL_2LVLS_Bind          ) = -1
    icntl( ICNTL_2LVLS_NNodes        ) = 0
    icntl( ICNTL_2LVLS_NcoresPerNode ) = 0
    icntl( ICNTL_2LVLS_NThrdsPerProc ) = 0
    icntl( ICNTL_2LVLS_NProcs        ) = 0

    ! Preconditioner
    icntl( ICNTL_PCD_Strategy ) = PCD_STRATEGY_isLocalApprox

    ! input system
    icntl( ICNTL_INSYSTEM) = INSYSTEM_isCentralized

    !---------------------------------------------------------------------------
    ! [2] Set dynamic default values
    !---------------------------------------------------------------------------

    ! select the first available sparse direct solver
    Do i = 1, SDS_MAX_INDEX

       Call XMPH_sds_select(i,sds_dummy,iinfo)

       if (iinfo >= 0) Then ! one found
          icntl(ICNTL_SDS_Default) = i
          Exit ! exit the loop
       End if

    End Do

    !---------------------------------------------------------------------------
    ! [3] End routine
    !---------------------------------------------------------------------------

    ! this routine do not return errors
    info = 0

  End Subroutine XMPH_maphys_Set_default_icntl


  ! [+] routine : XMPH_maphys_Set_default_rcntl  -------------------------------
  !
  !> initialize controls (reals)
  !!
  !! give default values to RCNTLs
  !!
  !! @param[in,out ] icntl   real controls
  !! @param[   out ] info    routine status
  !!
  !! @author Yohan Lee-tin-yien
  !! @todo   Implement
  !!
  Subroutine XMPH_maphys_Set_default_rcntl(rcntl,info)
    Use MPH_maphys_enum
    Implicit None

    !* Subroutine arguments *!

    real(kind=8), intent(inout) :: rcntl( MAPHYS_RCNTL_SIZE )
    integer, intent(  out) :: info

    !- End of header -----------------------------------------------------------

    ! Preconditioner
    rcntl(RCNTL_PCD_SparsifyThreshold) = 1.e-4

    ! Iterative solver 
    rcntl(RCNTL_ITS_Tolerance        ) = 1.e-5

    ! PILUT

    rcntl( RCNTL_ILU_LUThreshold    ) = 0.d0
    rcntl( RCNTL_ILU_SCHURThreshold ) = 0.d0

    info = 0

  End Subroutine XMPH_maphys_Set_default_rcntl


    ! [+] routine : XMPH_maphys_PretreatInputMatrix ---------------------------------
    !
    !> Pretreat the input matrix
    !!
    !!  Pretreat the matrix given by the user.
    !!  
    !!
    !!  @param[in,out] mphs   the maphys instance
    !!  @param[   in ] master the MPI rank where gb_A is defined
    !!  @param[   in ] rank   the MPI rank
    !!  @param[   out] gb_A   the pretreated input matrix 
    !!                        only relevant on "master"
    !!  @param[   out] info the routine status
    !!
    !!  @author Yohan Lee-tin-yien
    !!
    Subroutine XMPH_maphys_PretreatInputMatrix &
         (mphs, master, rank, gb_A, info )

      !* Modules *!
      Use MPH_maphys_enum
      Use MPH_error_mod
      Use MPH_mem_mod
      Use XMPH_maphys_type
      Use XMPH_sparse_matrix_mod
      Implicit None
      Include "mpif.h"

      !* Arguments *!
      Type(XMPH_maphys_t)       , Intent(inout) :: mphs
      Integer              , Intent(  in) :: master
      Integer              , Intent(  in) :: rank
      Type(XMPH_sparse_matrix_t), Intent(  out) :: gb_A
      Integer              , Intent(  out) :: info

      !* Local variables *!
      Real(kind=8) :: time_pretreatInputMatrix
      Real(kind=8) :: time_copyAssemble
      Real(kind=8) :: time_symStruct

      !- End of header ---------------------------------------------------------

      ! [-] Init

      info = 0
      mphs%rinfo (RINFO_TIMING_PrcInMat) = 0.d0 
      Call XMPH_sm_nullify(gb_A, info)
      MPH_ONFAILURE_RETURN(info)

      ! [--] Exit early
      If (rank /= master) Return

      time_pretreatInputMatrix = MPI_Wtime()

      ! [-] Copy inputs into gb_A
      !Write(*,*) "Warning:", Trim(FLNAME), __LINE__
      Call XMPH_sm_CreateFromData(&
           gb_A,SM_FMT_IJV,mphs%sym,&
           mphs%n, mphs%n, mphs%nnz, &
           mphs%rows, mphs%cols, mphs%values,&
           info )
      MPH_ONFAILURE_RETURN(info)
      
      ! [-] On symmetric matrices transpose the upper triangle
      !Write(*,*) "Warning:", Trim(FLNAME), __LINE__
      If ( (gb_A%sym == SM_SYM_IsSPD) .or. &
           (gb_A%sym == SM_SYM_IsSymmetric))Then
         Call XMPH_sm_transposeUpper(gb_A)
      End If
      
      ! [-] Assemble matrix
      !Write(*,*) "Warning:", Trim(FLNAME), __LINE__
      time_copyAssemble = MPI_Wtime()
      Call XMPH_sm_Assemble(gb_A, gb_A%n, info)
      CHCKASSRT( info >= 0, info )
      If (info < 0) Return
      time_copyAssemble = MPI_Wtime() - time_copyAssemble
      
      ! [-] Symetrize matrix structure on General matrices
      !Write(*,*) "Warning:", Trim(FLNAME), __LINE__
      If ( gb_A%sym == SM_SYM_IsGeneral )Then
         
         time_symStruct = MPI_Wtime()
         Call XMPH_sm_symStruct(gb_A, info)
         CHCKASSRT( info >= 0, info )
         If (info < 0) Return
         time_symStruct = MPI_Wtime() - time_symStruct
         
      End If
      
      ! [-] Convert the matrix into IJV format
      !Write(*,*) "Warning:", Trim(FLNAME), __LINE__
      Call XMPH_sm_convert( gb_A, SM_FMT_IJV, info )
      CHCKASSRT( info >= 0, info )
      If (info < 0) Return

      ! [-] write preprocessed matrix
      !Write(*,*) "Warning:", Trim(FLNAME), __LINE__
      If ( Trim(mphs%write_matrix) /= "" )Then
         
         Open(UNIT=111,FILE=mphs%write_matrix,&
              ACTION="WRITE",STATUS="REPLACE", IOSTAT=info)
         CHCKASSRT( info >= 0, info )
         If (info < 0) Return
         
         Call XMPH_sm_mmwrite(gb_A, 111, info )
         CHCKASSRT( info >= 0, info )
         If (info < 0) Return
         
      End If

      time_pretreatInputMatrix = MPI_Wtime() - time_pretreatInputMatrix

      ! Save statistics

      mphs%rinfo (RINFO_TIMING_PrcInMat) = time_pretreatInputMatrix
      mphs%iinfog(IINFOG_MAT_ORDER) = gb_A%m
      mphs%iinfog(IINFOG_MAT_NBENTRIES) = gb_A%nnz
      Call MPH_mem_add2usage(mphs%mem,XMPH_sm_sizeof( gb_A ))

    End Subroutine XMPH_maphys_PretreatInputMatrix

  End Module XMPH_maphys_aux_mod

