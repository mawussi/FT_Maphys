! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"

!> module for the linear system on the interface.
!!
!! 
Module XMPH_schur_mod

  !* Modules *!
  Use mph_error_mod
  Use XMPH_schur_assemble_mod, Only : &
       XMPH_schur_assemble_sparseMatrix
  Implicit None

  !* Private constants *!
  Character(len=MAPHYS_STRL), Private, Parameter :: &
       FLNAME = "xmph_schur_mod.F90"

  !* List of routines *!
  Public :: XMPH_SCHUR_InitTimers
  Public :: XMPH_SCHUR_Assemble_denseMatrix
  Public :: XMPH_SCHUR_Generate_rhs
  Public :: XMPH_schur_generate_initguess
  Public :: XMPH_SCHUR_GMRES_solve
  Public :: XMPH_SCHUR_FGMRES_solve 
  Public :: XMPH_SCHUR_CG_solve

Contains

  ! [+] routine : XMPH_SCHUR_InitTimers ----------------------------------
  !
  !> Initialize the timers of the iterative kernels 
  !!
  !! @param [in,out] mphs structure containing the timers
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_SCHUR_InitTimers(mphs)

    !* Modules *! 

    Use XMPH_maphys_type
    Use MPH_maphys_enum
    Implicit None

    !* Arguments *!

    Type(XMPH_maphys_t)       , Intent(inout) :: mphs

    !- End of header---------------------------------------------------------

    mphs%rinfo(RINFO_TIMING_SlvITSComm) = 0.d0
    mphs%rinfo(RINFO_TIMING_SlvITSMatV) = 0.d0
    mphs%rinfo(RINFO_TIMING_SlvITSPcdV) = 0.d0
    mphs%rinfo(RINFO_TIMING_Interp) = 0.d0
    mphs%rinfo(RINFO_RATIO_Lost) = 0.d0
    mphs%rinfo(RINFO_RATIO_Interp) = 0.d0
    mphs%rinfo(RINFO_PROC_INVOLVED) = 0.d0
5    mphs%rinfo(RINFO_TIMING_SlvITSDotP) = 0.d0

  End Subroutine XMPH_SCHUR_InitTimers

  ! [+] routine : XMPH_SCHUR_assemble_denseMatrix ------------------------------
  !
  !> Assemble a distributed dense Matrix defined on the interface.
  !!
  !! If a dense matrix is defined on the interface and distributed, 
  !! across multiple domain, this routine do :
  !!  
  !!       lc_MAT(i,j) = SUM lc_MAT(i,j) 
  !!
  !!
  !! @param[in     ] comm          MPI communicator
  !! @param[in     ] lc_domain     Specifies the local domain 
  !!
  !! @param[in,out ] lc_MAT        the local Matrix to assemble
  !!
  !! @param[out    ] info          the routine status
  !! 
  !!
  !! @author Azzam Haidar
  !! @author Luc   Giraud
  !! @author Yohan Lee-tin-yien
  !!
  !! @version 0.2
  !!
  !! @note
  !! From version 0.2, MPI_Bsend calls was substituted with MPI_Isend ones.
  !! 
  Subroutine XMPH_SCHUR_assemble_denseMatrix( & ! intents
       comm, lc_domain,                   & ! in
       lc_MAT,                            & ! inout
       info                               & ! out
       )

    !* Module(s) *!
    Use XMPH_maphys_type, Only : &
         maphys_domain_t, & ! types
         XMPH_dense_matrix_t
    Use mph_domain_mod
    Implicit None
    Include 'mpif.h'

    !* Arguments *!
    Integer               , Intent(in   ) :: comm        
    Type(maphys_domain_t) , Intent(in   ) :: lc_domain        
    Type(XMPH_dense_matrix_t)  , Intent(inout) :: lc_MAT    
    Integer               , Intent(out  ) :: info         

    !* Local variables *!

    ! Scalars
    Integer, Parameter :: PresetMPITagg = 77
    Integer :: tagg, rank 
    Integer :: sizeIntrf
    Integer :: nNeighb 
    Integer :: ldA
    Integer :: neighb  
    MPH_INT :: maxSizeIntrf
    MPH_INT :: recvbuffSize, sendBuffSize
    MPH_INT :: sendBuffStart
    MPH_INT :: i,j,ij,k,l          

    ! Arrays 
    Integer     , Pointer :: indexVi          (:)
    Integer     , Pointer :: ptr_Index_Intrfc (:)
    Integer     , Pointer :: Index_Intrfc     (:)
    XMPH_FLOAT, Pointer :: A                (:)
    XMPH_FLOAT, Pointer :: sendbuff             (:)
    XMPH_FLOAT, Pointer :: recvbuff             (:)

    Integer :: mpistatus(MPI_STATUS_SIZE) 
    Integer, Pointer :: mpireq(:) ! MPI Request

    !- End of header ---------------------------------------------------------

    !azz> begin previous header 
    ! s                      : local Schur matrix of this subdomain to be completed
    !                          with its neighbors contributions.
    !
    ! sizeIntrf              : size of the interface (and of the Schur matrix)
    !
    ! nNeighb                   : number of subdomains neighbor of the subdomain
    !                          no (me+1) as me is the MPI process number starting
    !                          at 0.
    !
    ! indexVi(nNeighb)          : no of the neighboring subdomain
    !
    ! ptr_Index_Intrfc(nNeighb) : pointer on the index of the first point on the
    !                          interface with the sudomain. The indices of the
    !                          interface points are stored in the array 
    !                          Index_Intrfc.
    !                          The indices are given in the local ordering.
    !
    ! Index_Intrfc(*)        : list of the interface points, pointed by the
    !                          array ptr_Index_Intrfc.
    !                          The indices are given in the local ordering.
    !
    !azz< end previous header   

    !-------------------------------------------------------------------------
    ! [1] Initialize local variables
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! [1.1] Nullify all arrays
    !-------------------------------------------------------------------

    ! Nullify All arrays

    Nullify( indexVi, ptr_Index_Intrfc, Index_Intrfc )
    Nullify( A )
    Nullify( sendbuff, recvbuff, mpiReq )

    !-------------------------------------------------------------------
    ! [1.2] Check arguments
    !-------------------------------------------------------------------

    ! comm

    info = 0
    rank  = 0
    tagg  = PresetMPITagg

    Call MPI_Comm_rank(comm,rank,info)
    ASSRT( info == MPI_SUCCESS )

    ! lc_domain
    Call MPH_domain_check( lc_domain, info )
    CHCKASSRT( info >= 0, info )
    If ( info < 0 ) Goto 9999 

    ! lc_MAT
    CHCKASSRT( Associated( lc_MAT%v ), info )
    If( info < 0 ) Goto 9999

    !-------------------------------------------------------------------
    ! [1.3] Real values
    !-------------------------------------------------------------------

    ! Domain related

    sizeIntrf = lc_domain%mysizeIntrf
    nNeighb      = lc_domain%mynbvi

    indexVi          => lc_domain%myindexVi
    ptr_Index_Intrfc => lc_domain%myptrindexVi 
    Index_Intrfc     => lc_domain%myindexIntrf

    ! matrix related

    A   => lc_MAT%v
    ldA =  lc_MAT%ld

    !-------------------------------------------------------------------------
    ! [2] Allocate the communication buffer
    !-------------------------------------------------------------------------

    ! Get the sizes 

    maxSizeIntrf = 0
    sendBuffSize = 0

    Do neighb = 1,nNeighb

       ! k is the size of the interface, 
       ! sqr(k) the size of the portion of the Schur we send
       k = ptr_Index_Intrfc(neighb+1)-ptr_Index_Intrfc(neighb)
       if ( k > maxSizeIntrf ) maxSizeIntrf = k
       sendBuffSize = sendBuffSize + k**2

    End Do
    recvbuffSize =  maxSizeIntrf**2 

    ! Perform allocations

    Allocate ( &
         recvbuff(recvbuffSize), &
         sendBuff(sendBuffSize), &
         mpiReq(nNeighb), &
         STAT = info )

    CHCKASSRT( info == 0, info )
    If( info < 0 ) Goto 9999

    !-------------------------------------------------------------------------
    ! [3] Send this subdomain's contributions to each neighbor's Schur matrix 
    !-------------------------------------------------------------------------

    sendBuffStart = 1
    Do neighb =1,nNeighb

       ! Copy values to the buffer
       k = sendBuffStart
       Do j=ptr_Index_Intrfc(neighb), ptr_Index_Intrfc(neighb+1)-1
          Do i=ptr_Index_Intrfc(neighb), ptr_Index_Intrfc(neighb+1)-1

             ij = Index_Intrfc(i) + ldA*( Index_Intrfc(j) - 1)
             sendbuff(k) = A(ij)
             k = k + 1

          End Do
       End Do

       ! k is the size of the interface, 
       ! sqr(k) the size of the portion of the Schur we send
       k = ptr_Index_Intrfc(neighb+1)-ptr_Index_Intrfc(neighb)
       k = k*k
       ! use a buffered send instead of a mpi_send because of problems with
       ! large messages 
       CALL MPI_Isend(                    &
            sendbuff(sendBuffStart),k,XMPH_FLOATMPI,indexVi(neighb),&
            tagg,comm,mpiReq(neighb),info)
       ASSRT( info == MPI_SUCCESS )

       sendBuffStart = sendBuffStart + k
    End Do

    !-------------------------------------------------------------------------
    ! [4] Receive the values from neighbors
    !-------------------------------------------------------------------------

    Do l=1,nNeighb

       ! Get a message
       Call MPI_Probe(MPI_ANY_SOURCE,tagg,comm,mpiStatus,info)
       ASSRT( info == MPI_SUCCESS )
       
       ! Detect which neighbor this message correspond
       
       !lty> if the number of neighbors is too big, 
       !lty> we should store the source2neigbIndex in a table.
       Do neighb=1,nNeighb
          If ( indexVi(neighb)  ==  mpiStatus(MPI_SOURCE) ) Exit
       End Do

       ! Receive the buffer
       k = ptr_Index_Intrfc(neighb+1)-ptr_Index_Intrfc(neighb)
       k = k*k
       Call MPI_Recv(recvbuff(1),k,&
            XMPH_FLOATMPI,mpiStatus(MPI_SOURCE),&
            tagg,comm,MPI_STATUS_IGNORE,info)
       ASSRT( info == MPI_SUCCESS )

       ! Add contribution 
       k = 1                     
       Do j=ptr_Index_Intrfc(neighb), ptr_Index_Intrfc(neighb+1) -1
          Do i=ptr_Index_Intrfc(neighb), ptr_Index_Intrfc(neighb+1) -1

             ij = Index_Intrfc(i) + ldA*( Index_Intrfc(j) - 1)
             A(ij) = A(ij) + recvbuff(k)
             k = k + 1

          End Do
       End Do

    End Do

    ! Wait for all MPI_Isend to complete
    Call MPI_Waitall(nNeighb, mpiReq, MPI_STATUSES_IGNORE,info)
    ASSRT( info == MPI_SUCCESS )
    info = 0

    !-------------------------------------------------------------------------
    ! [5] Exit routine
    !-------------------------------------------------------------------------
    
9999 Continue

    ! Free memory
    If (Associated(mpiReq)) Deallocate (mpiReq)
    If (Associated(sendbuff)) Deallocate (sendbuff)
    If (Associated(recvbuff)) Deallocate (recvbuff)


   End Subroutine XMPH_SCHUR_assemble_denseMatrix

  ! [+] routine : XMPH_SCHUR_Generate_InitGuess ----------------------------------
  ! 
  !> Generate the Initial Guess for the interface system.
  !!
  !! Generate the initial guess for the interface system "bound_sol" such as :
  !!
  !! bound_sol := SUM [ domain_SOLb ]
  !!
  !!***** 
  !! 
  !! @param[in,out] mphs 
  !!
  !!       The maphys instance containing :
  !!       - the necessary data to SUM on all the processors
  !!
  !!       This routine do not append any data to mphs. It only modify
  !!       temporary buffers and arrays. Hence, mphs is conceptualy an "intent(in)".
  !!
  !! @param[in,out] domain_sol
  !!
  !!       The init guess of the domain linear system.
  !!
  !! @param[in,out] bound_sol
  !!
  !!       The init guess of the interface linear system, distributed.
  !!
  !! @param[in,out] info
  !!
  !!       The routine status
  !!
  !!***** 
  !!
  !! @author Yohan Lee-tin-yien
  !!
  !!
  !! @see For details, see Azzam's Thesis, section 2.13, page 19, formula 2.7.
  !!
  Subroutine XMPH_SCHUR_Generate_InitGuess( & ! intents
       mphs,                           & ! inout
       domain_sol,                     & ! in
       bound_sol,                      & ! out
       info                            &
       )

    !* Modules *!

    Use XMPH_maphys_type
    Use XMPH_dense_matrix_mod
    Use XMPH_schur_aux_mod
    Use mph_error_mod
    Implicit None

    !* Arguments *!

    Type(XMPH_maphys_t      ), Intent(inout) :: mphs
    Type(XMPH_dense_matrix_t), Intent(in   ) :: domain_sol
    Type(XMPH_dense_matrix_t), Intent(  out) :: bound_sol
    Integer             , Intent(  out) :: info

    !* Local variables *!
    MPH_INT                     :: bound_size
    MPH_INT                     :: ioffset

    !----------------------------------------------------------------------
    !  [1] Allocate structure
    !----------------------------------------------------------------------
    
    bound_size = mphs%lc_domain%mysizeIntrf
    Call XMPH_dm_create &
         (bound_sol,bound_size,1,bound_size, info) 
    CHCKASSRT( info >= 0, info)
    If (info < 0 ) Return

    !----------------------------------------------------------------------
    !  [2] Get domain_solB (saved in bound_sol)
    !----------------------------------------------------------------------
    
    ioffset= mphs%lc_domain%myndof - bound_size
    Call XMPH_dm_copyBloc(                  &
         bound_sol%m  ,bound_sol%n,              &
         domain_sol%ld,&
         domain_sol%v(ioffset+1:ioffset+domain_sol%ld*bound_sol%m), &
         bound_sol%ld ,&
         bound_sol%v(1:bound_sol%ld*bound_sol%m), &
         info)
    CHCKASSRT( info >= 0, info)
    If (info < 0 ) Return

    !----------------------------------------------------------------------
    !  [3] Sum on all processors
    !----------------------------------------------------------------------
    
    Call XMPH_SCHUR_UpdateVector(mphs, bound_sol , info )
    CHCKASSRT( info >= 0, info)
    If (info < 0 ) Return

  End Subroutine XMPH_SCHUR_Generate_InitGuess





  ! [+] routine : XMPH_SCHUR_Generate_RHS -----------------------------------------
  ! 
  !> Generate the RHS for the interface system.
  !!
  !! Generate the RHS for the interface system "bound_rhs" such as :
  !!
  !! 1) bound_rhs := SUM [ domain_RHSb - (Abi . Aii^-1 . domain_RHSi ) ]
  !! 
  !!    or 
  !!
  !! 2) bound_rhs := domain_RHSb - SUM [ (Abi . Aii^-1 . domain_RHSi ) ]
  !!
  !! depending on how the global RHS was distributed.
  !!
  !! Currently, the formula 2) is forced.
  !! 
  !!***** 
  !! 
  !! @param[in,out] mphs 
  !!
  !!       The maphys instance containing :
  !!       - the factors of Aii
  !!       - the matrix     Abi
  !!       - the necessary data to SUM on all the processors
  !!
  !!       This routine do not append any data to mphs. It only modify
  !!       temporary buffers and arrays. Hence, mphs is conceptualy an "intent(in)".
  !!
  !! @param[in,out] domain_rhs
  !!
  !!       The RHS of the domain linear system.
  !!
  !! @param[in,out] bound_rhs
  !!
  !!       The local part of the RHS of the interface linear system, distributed.
  !!
  !! @param[in,out] info
  !!
  !!       The routine status
  !!
  !!***** 
  !!
  !! @version 0.1    Yohan Lee-tin-yien
  !!        
  !!     Refactorize the routine :
  !!     modify input/output, variables names, etc.
  !!
  !! @version alpha  Azzam 
  !!
  !!     initial implementation
  !!
  Subroutine XMPH_SCHUR_Generate_RHS( & ! intents
       mphs,                      & ! inout
       domain_rhs,                & ! in
       bound_rhs,                & ! out
       info                       &
       )

    !* Modules *!

    Use XMPH_maphys_type
    Use XMPH_sds_mod
    Use XMPH_dense_matrix_mod
    Use XMPH_schur_aux_mod
    Use mph_log_mod
    Use MPH_maphys_enum, Only : &
         IINFOG_ITS_NbIters        ! enumerations
    Implicit None

    !* Arguments *!

    Type(XMPH_maphys_t      ), Intent(inout) :: mphs
    Type(XMPH_dense_matrix_t), Intent(in   ) :: domain_rhs
    Type(XMPH_dense_matrix_t), Intent(  out) :: bound_rhs
    Integer             , Intent(  out) :: info

    !* Local variables *!

    ! constants
    Integer, Parameter :: rhsway = 2

    ! Scalars
    MPH_INT :: i
    MPH_INT :: lim
    MPH_INT :: domain_ndof
    MPH_INT :: bound_ndof

    ! Derived types
    Type(XMPH_dense_matrix_t) :: temp

    !----------------------------------------------------------------------
    !  [1] Initialize variables
    !----------------------------------------------------------------------

    ! scalars
    ! rhsway      = mphs%part_rhs%rhsway ! lty: forced option
    domain_ndof = mphs%lc_domain%myndof
    bound_ndof  = mphs%lc_domain%mysizeIntrf
    lim         = domain_ndof - bound_ndof

    ! result
    Call XMPH_dm_create(bound_rhs,bound_ndof,1,bound_ndof,info)   
    MPH_ONFAILURE_GOTO9999(info)

    !----------------------------------------------------------------------
    !  [2] compute : temp <- Aii^{-1} . domain_RHSi
    !----------------------------------------------------------------------

    ! Setup the right hand side of the Schur system

    Call XMPH_dm_dup(domain_rhs,temp, info )   
    MPH_ONFAILURE_GOTO9999(info)

    temp%v(lim+1:domain_ndof) = XMPH_FLOATZERO


    ! Compute Aii^{-1}.domain_RHSi

    Call XMPH_sds_solve_RHS(             &
         mphs%sls_domain%sds,       & ! temp = AII^{-1}.bI
         temp%v, temp%n, temp%ld, &
         info )
    MPH_ONFAILURE_GOTO9999(info)
    
    !----------------------------------------------------------------------
    !  [3] compute : bound_rhs <-- Abi . temp
    !----------------------------------------------------------------------

    Call XMPH_sm_vectorproduct(mphs%sm_Abi, -lim, 0,temp, bound_rhs, info )
    MPH_ONFAILURE_GOTO9999(info)

    !----------------------------------------------------------------------
    !  [4] compute : bound_rhs <-- SUM[ domain_RHSb - bound_rhs ]
    !  [ - ] or        bound_rhs <-- domain_RHSb - SUM[ bound_rhs ]
    !  [ - ] depending on strategy to distribute the RHS
    !----------------------------------------------------------------------

    !----------------------------------------------------------------------
    !  [4.1] RHSWAY == 1 : compute bound_rhs <-- SUM( bound_rhs )
    !----------------------------------------------------------------------

    !azz> 
    !     If rhs is already update or assembled so just 
    !     assemble the contribution Abi * Aii^{-1}* bI then subtract from bB 
    !     otherwise do bB- ABI * AII^{-1}* bI then assemble the result 
    !    
    !azz<

    If( RHSWAY == 1 )Then 

       Call XMPH_SCHUR_UpdateVector(mphs, bound_rhs , info )
       MPH_ONFAILURE_GOTO9999(info)

    Endif

    !----------------------------------------------------------------------
    !  [4.2] RHSWAY == 1 : compute bound_rhs <-- domain_RHSb - bound_rhs
    !----------------------------------------------------------------------

    Do i= 1, bound_rhs%m
       bound_rhs%v(i) = domain_rhs%v(lim+i) - bound_rhs%v(i)
    End Do

    !----------------------------------------------------------------------
    !  [4.3] RHSWAY == 2 : compute  bound_rhs <-- SUM( bound_rhs )
    !----------------------------------------------------------------------
    !azz>
    !     Assemble the right hand side on the boundaries
    !     IF RHS IS ALREADY UPDATE OR ASSEMBLED SO JUST 
    !     ASSEMBLE THE CONTRIBUTION ABI * AII^{-1}* bI THEN SUBTRACT FROM bB 
    !     OTHERWISE DO bB- ABI * AII^{-1}* bI THEN ASSEMBLE THE RESULT 
    !azz<

    If( RHSWAY == 2 )Then 

       Call XMPH_SCHUR_UpdateVector(mphs, bound_rhs , info )
       MPH_ONFAILURE_GOTO9999(info)

    Endif

    !----------------------------------------------------------------------
    !  [5] Exit routine
    !----------------------------------------------------------------------

9999 Continue

    ! Free memory
    Call XMPH_dm_free(temp,info)

  End Subroutine XMPH_SCHUR_Generate_RHS




  ! [+] routine :  XMPH_SCHUR_GMRES_solve -----------------------------------------
  !
  !> solve the interface linear system using GMRES
  !!
  !! Solve the Schur complement system (distributed) with GMRES.
  !! where the system is :
  !!           S. x = rhs
  !!
  !!-----
  !!
  !! @param[in,out] mphs  
  !!
  !!       the maphys instance which contains :
  !!       - the parameters of the Iterative solver
  !!       - the Schur "S"      (which have multiple representation)
  !!       - Its Preconditioner (which have multiple representation)
  !!       - necessary data to communicate between processors.
  !!       - the informations : iinfo/rinfo [in,out]
  !!
  !! @param[in    ] bound_rhs  
  !!
  !!       it holds the local part of the right-hand-side "rhs"
  !!
  !! @param[in,out] bound_sol
  !!
  !!       - On input , it may hold the local part of the initial guess x0
  !!       - On output, it holds the local part of the solution x
  !!
  !! @param[   out] info 
  !!       
  !!       Integer describing the routine status
  !!
  !!-----
  !!
  !! @note
  !! See also the documentation of PackGMRES at :
  !!   http://www.cerfacs.fr/algor/reports/2003/TR_PA_03_03.pdf
  !!
  !! @author Azzam Haidar
  !! @author Luc   Giraud
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_SCHUR_GMRES_solve(mphs, bound_rhs, bound_sol, info)

    !* Module used *!

    Use XMPH_maphys_type
    Use MPH_maphys_enum
    Use XMPH_schur_aux_mod
    Use XMPH_dense_matrix_mod
    Use XMPH_fault_mod
    Use mph_log_mod
    Use mph_error_mod
    Use MPH_maphys_enum, Only : &
         IINFOG_ITS_NbIters        ! enumerations
    Use MPH_mem_mod
    Implicit None

    !* Arguments *!

    Type(XMPH_maphys_t)        , Intent(inout) :: mphs
    Type(XMPH_dense_matrix_t)  , Intent(in   ) :: bound_rhs 
    Type(XMPH_dense_matrix_t)  , Intent(inout) :: bound_sol 
    Integer               , Intent(  out) :: info

    !* Local Variables *!

    ! Scalars
    Integer    :: lworkStrategy
    Integer    :: restrt
    MPH_INT :: lwork
    MPH_INT :: intrf_ndof
    MPH_INT :: bound_ndof
    MPH_INT :: interior_ndof
    MPH_INT :: domain_ndof
    Integer    :: i
    Integer    :: fault_info
    Integer    :: rank, master
    Integer    :: verb
    Logical    :: HasInitGuess
    Integer(kind=8) :: mem

    ! Arrays
    Integer, Target :: gmirc   (5)
    Integer         :: gmicntl (8)
    Real(Kind=XMPH_FLOATKIND)         :: gmrcntl (5)

    Integer         :: gmiinfo (8)
    Real(Kind=XMPH_FLOATKIND)         :: gmrinfo (4)    !lty: why 4 and not 1 ?

    XMPH_FLOAT, Pointer :: work     (:)

    ! Aliases (Pointers)

    ! to gmirc (reverse communication)
    Integer, Pointer :: revcom
    Integer, Pointer :: colx
    Integer, Pointer :: coly
    Integer, Pointer :: colz
    Integer, Pointer :: nbscal

    ! to gmres workspace "work"
    Type(XMPH_dense_matrix_t) :: x 
    Type(XMPH_dense_matrix_t) :: y 
    Type(XMPH_dense_matrix_t) :: z

    !- End of header -------------------------------------------------

    !-----------------------------------------------------------------
    ! [1] Initialize local variables / GMRES package
    !-----------------------------------------------------------------

    !---------------------------------------------------------
    ! [1.1] Get differents sizes & options
    !---------------------------------------------------------

    intrf_ndof    = mphs%lc_domain%gballintrf
    bound_ndof    = mphs%lc_domain%myndofintrf
    interior_ndof = mphs%lc_domain%myndofinterior
    domain_ndof   = mphs%lc_domain%myndof
    HasInitGuess = (mphs%IKEEP(IKEEP_ITS_InitGuess) == 1)

    !---------------------------------------------------------
    ! [1.2] Set GMRES options
    !---------------------------------------------------------


    ! Get default options
    Call init_XMPH_ARITHgmres(gmicntl,gmrcntl)

    ! Set specific options
    gmrcntl(1) = Real( mphs%rkeep(RKEEP_ITS_Tolerance),KIND=XMPH_FLOATKIND)

    gmicntl(2) = 0
    gmicntl(5) = mphs%ikeep(IKEEP_ITS_OrtStrategy )
    gmicntl(6) = mphs%ikeep(IKEEP_ITS_InitGuess   )
    gmicntl(7) = mphs%ikeep(IKEEP_ITS_MaxIters    )
    gmicntl(8) = mphs%ikeep(IKEEP_ITS_ResStrategy )

    restrt     = mphs%ikeep(IKEEP_ITS_Restart     )

    ! Outputs

    rank= mphs%ikeep(IKEEP_MPIRANK)
    master= mphs%ikeep(IKEEP_HOSTRANK)
    Call mph_log_GetVerbosity(verb)
    If ( (verb >= MSG_STD).and.(rank == master) )Then
       ! print warning messages
       Call mph_log_GetUnit(MSG_WARNING, gmicntl(2))
       ! print convergence history to stdout
       Call mph_log_GetUnit(MSG_STD    , gmicntl(3))
    End If


    ! Force right preconditioning we use a preconditioner
    If ( mphs%ikeep(IKEEP_Pcd_Strategy ) == PCD_STRATEGY_isNone ) &
         gmicntl(4) = 0
    If ( mphs%ikeep(IKEEP_Pcd_Strategy ) /= PCD_STRATEGY_isNone ) &
         gmicntl(4) = 2

    ! Create solution if no init guess
    If (HasInitGuess .eqv. .False.) Then
       Call XMPH_dm_Nullify(bound_sol,info) 
       MPH_ONFAILURE_GOTO9999(info)
       Call XMPH_dm_Create(bound_sol,bound_ndof,1,bound_ndof,info  ) 
       MPH_ONFAILURE_GOTO9999(info)
       Do i=1,bound_sol%m
          bound_sol%v(i) = XMPH_FLOATZERO
       End Do
    End If

    !---------------------------------------------------------
    ! [1.3] Allocate GMRES workspace "work"
    !---------------------------------------------------------

    ! select which formula to compute lwork
    If (( gmicntl(5) == 0) .or. (gmicntl(5) == 1)) Then
       If ( gmicntl(8) == 1) lworkStrategy = 1
       If ( gmicntl(8) /= 1) lworkStrategy = 2
    Else
       If ( gmicntl(8) == 1) lworkStrategy = 3
       If ( gmicntl(8) /= 1) lworkStrategy = 4
    Endif

    ! Set lwork value
    lwork = restrt**2 + restrt*(bound_ndof+5) + 5* bound_ndof + 1
    Select Case(lworkStrategy)
    Case (1); lwork = lwork + (         1        )    
    Case (2); lwork = lwork + (bound_ndof + 1     )
    Case (3); lwork = lwork + (restrt            )
    Case (4); lwork = lwork + (bound_ndof + restrt)
    End Select

    ! Allocate
    
    Allocate(work(lwork), STAT = info )
    CHCKALLOC(info)
    MPH_ONFAILURE_GOTO9999(info)
    
    !---------------------------------------------------------
    ! [1.4] Setup the data
    !---------------------------------------------------------


    ! Setup the rhs
    ! Work(1:bound_ndof)             = local initial guess
    ! Work(bound_ndof+1:2*bound_ndof) = local right-hand-side
    If (HasInitGuess) Call XMPH_ARITHcopy &
         (bound_ndof, bound_sol%v(1),1, work(1), 1)

    Call XMPH_ARITHcopy&
         (bound_ndof,bound_rhs%v(1),1,work(bound_ndof+1), 1)

    ! initialize vectors
    info = 0
    Call XMPH_dm_nullify(x,info)
    Call XMPH_dm_nullify(y,info)
    Call XMPH_dm_nullify(z,info)
    MPH_ONFAILURE_GOTO9999(info)

    x%ld = bound_ndof
    x%m  = bound_ndof
    x%n  = 1

    y%ld = bound_ndof
    y%m  = bound_ndof
    y%n  = 1

    z%ld = bound_ndof
    z%m  = bound_ndof
    z%n  = 1

    ! Associate aliases
    revcom => gmirc(1)
    colx   => gmirc(2)
    coly   => gmirc(3)
    colz   => gmirc(4)
    nbscal => gmirc(5)

    !-----------------------------------------------------------------
    ! [2] Call the solver, Iterate until solution is found
    !-----------------------------------------------------------------

    revcom = 1 
    Do While ((revcom > 0).And.(info >= 0))

       !----------------------------------------------------
       ! [2.1] Call the driver
       !----------------------------------------------------
       gmiinfo(2) = 0  !fault 
       gmirc(8) =0       ! fault
       Call Drive_XMPH_ARITHgmres(         &
            intrf_ndof,bound_ndof,  &
            restrt,lwork,work,    &
            gmirc,gmicntl,gmrcntl, &
            gmiinfo,gmrinfo)

       !----------------------------------------------------
       ! [2.2] Associate arrays
       !----------------------------------------------------

       y%v => work(coly : (coly+bound_ndof-1) )
       z%v => work(colz : (colz+bound_ndof-1) )

       If ( revcom /= 4)Then
          x%n   = 1
          x%v => work(colx : (colx+bound_ndof-1) )
       Else
          x%n   =  nbscal
          x%v => work(colx : (colx + bound_ndof*nbscal -1) )
       End If

       !----------------------------------------------------
       ! [2.2.1] Fault simulation
       !----------------------------------------------------
!        !check if the iteration will be faulty
!         If(gmirc(8) .Eq. 0) Then
!           Call XMPH_IsFault(mphs, gmiinfo(2),fault_info)
!           If(fault_info .Eq. 2) Then !multiple faults not handle in GMRES
!              gmirc(8) = 1
!           Else If(fault_info .Eq. 1)Then !single fault
!              Call XMPH_Check_Fault(mphs,lwork, work, gmiinfo(2),gmicntl)              
!           End If
!        End If

       !----------------------------------------------------
       ! [2.3] Do what the driver asked to perform.
       !----------------------------------------------------
       !
       ! 1     - Perform Matrix Vector Product : z <-- Schur . x
       ! 2     - Apply left preconditioning    : z <-- Precond . x
       ! 3     - Apply right precondiotioning  : z <-- Precond . x
       ! 4     - Perform scalar product        : z <-- x * y
       ! 0     - Exit                          : Solving step ended
       ! other - Throw Error : unknown revcom
       !
       
       Select Case (revcom)
       Case (0); Continue                      ! Exit    
       Case (1);  Call XMPH_SCHUR_SchurVectorProduct (mphs, x, z, info)
       Case (2); CHCKASSRT(.False., info )     ! MaPHyS does not support left preconditioning
       Case (3); Call XMPH_SCHUR_PcdVectorProduct   (mphs, x,    z, info)
       Case (4); Call XMPH_SCHUR_VectorScalarProduct(mphs, x, y, z, info)
       Case Default; CHCKASSRT(.False., info ) ! Undocumented revcom
       End Select

    End Do

    !---------------------------------------------------------------------------
    ! [-] Checking the outputs 
    !---------------------------------------------------------------------------    

    !-----------------------------------------------
    ! [--] Checking revcom, print a message on error
    !-----------------------------------------------

    Select Case (revcom)
    Case (0); Continue ! Gmres exited correctly
    Case (1); Call MPH_Log(MSG_ERROR,"Error occured in matrix/vector Product")
    Case (2); Call MPH_Log(MSG_ERROR,"MaPHyS does not support left preconditioning")
    Case (3); Call MPH_Log(MSG_ERROR,"Error occured in right preconditioning")
    Case (4); Call MPH_Log(MSG_ERROR,"Error occured in scalar product")
    Case Default; Call MPH_LogWithInfo(MSG_ERROR,revcom,&
            "From PackGMRES: undocumented revcom code. revcom =")
    End Select
    MPH_ONFAILURE_GOTO9999(info)
    
    !----------------------------------------------------
    ! [--] Checking gmiinfo(1), the return status of PackGMRES
    !----------------------------------------------------

    Select Case (gmiinfo(1))
    Case ( 0) ! convergence achieved
       info = MPH_SUCCESS
       Continue
    Case (-4) ! non convergence. Consider it as a warning and not an error
       info = MPH_SUCCESS + 1
       Call MPH_LogWithInfo(MSG_WARNING,gmicntl(7),&
            "From PackGMRES: convergence not achieved. Nb iterations:")
    Case Default
       Call MPH_LogWithInfo(MSG_ERROR,gmiinfo(1),"PackGRMES exited with error code =")
       CHCKASSRT(.False., info )
    End Select
    MPH_ONFAILURE_GOTO9999(info)

    !-----------------------------------------------------------------
    ! [-] Exit routine
    !-----------------------------------------------------------------

    ! Save data
    Call XMPH_ARITHcopy(bound_ndof, work(1),1, bound_sol%v(1),1 )

    mphs%IINFOG( IINFOG_ITS_NbIters ) = gmiinfo(2)
    mphs%RINFOG( RINFOG_ITS_BckErr  ) = gmrinfo(1)

    mem = INT(lwork*XMPH_FLOATBYTESIZE,8)
    Call MPH_mem_add2usage(mphs%mem, mem )
    mphs%IINFO ( IINFO_ITS_MEMPEAK  ) = byte2Mbyte(mem)
    mphs%IINFO ( IINFO_ITS_MEMUSED  ) = byte2Mbyte(mem)

    ! Free memory
9999 Continue
    If (Associated(work)) Deallocate( work )

  End Subroutine XMPH_SCHUR_GMRES_Solve



  ! [+] routine :  XMPH_SCHUR_FGMRES_solve -----------------------------------------
  !
  !> solve the interface linear system using FGMRES
  !!
  !! Solve the Schur complement system (distributed) with FGMRES.
  !! where the system is :
  !!           S. x = rhs
  !!
  !!-----
  !!
  !! @param[in,out] mphs  
  !!
  !!       the maphys instance which contains :
  !!       - the parameters of the Iterative solver
  !!       - the Schur "S"      (which have multiple representation)
  !!       - Its Preconditioner (which have multiple representation)
  !!       - necessary data to communicate between processors.
  !!       - the informations : iinfo/rinfo [in,out]
  !!
  !! @param[in    ] bound_rhs  
  !!
  !!       it holds the local part of the right-hand-side "rhs"
  !!
  !! @param[in,out] bound_sol
  !!
  !!       - On input , it may hold the local part of the initial guess x0
  !!       - On output, it holds the local part of the solution x
  !!
  !! @param[   out] info 
  !!       
  !!       Integer describing the routine status
  !!
  !!-----
  !!
  !! @note
  !! See also the documentation of PackFGMRES at :
  !!   http://www.cerfacs.fr/algor/reports/2003/TR_PA_03_03.pdf
  !!
  !! @author Azzam Haidar
  !! @author Luc   Giraud
  !! @author Yohan Lee-tin-yien
  !! @author Mawussi Zounon 
  Subroutine XMPH_SCHUR_FGMRES_solve(mphs, bound_rhs, bound_sol, info)

    !* Module used *!

    Use XMPH_maphys_type
    Use MPH_maphys_enum
    Use XMPH_schur_aux_mod
    Use XMPH_dense_matrix_mod
    Use XMPH_fault_mod
    Use mph_log_mod
    Use mph_error_mod
!    Use MPH_maphys_enum, Only : &
!         IINFOG_ITS_NbIters        ! enumerations
    Use MPH_mem_mod
    Implicit None

    !* Arguments *!

    Type(XMPH_maphys_t)        , Intent(inout) :: mphs
    Type(XMPH_dense_matrix_t)  , Intent(in   ) :: bound_rhs 
    Type(XMPH_dense_matrix_t)  , Intent(inout) :: bound_sol 
    Integer               , Intent(  out) :: info

    !* Local Variables *!

    ! Scalars
    Integer    :: restrt
    MPH_INT :: lwork
    MPH_INT :: intrf_ndof
    MPH_INT :: bound_ndof
    MPH_INT :: interior_ndof
    MPH_INT :: domain_ndof
    Integer    :: i
    Integer    :: fault_info, last_failed_it
    Integer    :: rank, master
    Integer    :: verb
    Logical    :: HasInitGuess
    Integer(kind=8) :: mem

    ! Arrays
    Integer, Target :: fgmirc   (8)
    Integer         :: fgmicntl (6)
    Real(Kind=XMPH_FLOATKIND)         :: fgmrcntl (3)

    Integer         :: fgmiinfo (8)
    Real(Kind=XMPH_FLOATKIND)         :: fgmrinfo     !lty: why 4 and not 1 ?
    
    XMPH_FLOAT, Pointer :: work     (:)

    ! Aliases (Pointers)

    ! to fgmirc (reverse communication)
    Integer, Pointer :: revcom
    Integer, Pointer :: colx
    Integer, Pointer :: coly
    Integer, Pointer :: colz
    Integer, Pointer :: nbscal

    ! to fgmres workspace "work"
    Type(XMPH_dense_matrix_t) :: x 
    Type(XMPH_dense_matrix_t) :: y 
    Type(XMPH_dense_matrix_t) :: z

    !- End of header -------------------------------------------------

    !-----------------------------------------------------------------
    ! [1] Initialize local variables / FGMRES package
    !-----------------------------------------------------------------

    !---------------------------------------------------------
    ! [1.1] Get differents sizes & options
    !---------------------------------------------------------

    intrf_ndof    = mphs%lc_domain%gballintrf
    bound_ndof    = mphs%lc_domain%myndofintrf
    interior_ndof = mphs%lc_domain%myndofinterior
    domain_ndof   = mphs%lc_domain%myndof
    HasInitGuess = (mphs%IKEEP(IKEEP_ITS_InitGuess) == 1)

    !---------------------------------------------------------
    ! [1.2] Set FGMRES options
    !---------------------------------------------------------


    ! Get default options
    Call init_XMPH_ARITHfgmres(fgmicntl,fgmrcntl)

    ! Set specific options
    fgmrcntl(1) = Real( mphs%rkeep(RKEEP_ITS_Tolerance),KIND=XMPH_FLOATKIND)

    fgmicntl(2) = 0
    fgmicntl(4) = mphs%ikeep(IKEEP_ITS_OrtStrategy )
    fgmicntl(5) = mphs%ikeep(IKEEP_ITS_InitGuess   )
    fgmicntl(6) = mphs%ikeep(IKEEP_ITS_MaxIters    )

    restrt     = mphs%ikeep(IKEEP_ITS_Restart     )

    ! Outputs

    rank= mphs%ikeep(IKEEP_MPIRANK)
    master= mphs%ikeep(IKEEP_HOSTRANK)
    Call mph_log_GetVerbosity(verb)
    If ( (verb >= MSG_STD).and.(rank == master) )Then
       ! print warning messages
       Call mph_log_GetUnit(MSG_WARNING, fgmicntl(2))
       ! print convergence history to stdout
       Call mph_log_GetUnit(MSG_STD    , fgmicntl(3))
    End If

    ! Create solution if no init guess
    If (HasInitGuess .eqv. .False.) Then
       Call XMPH_dm_Nullify(bound_sol,info) 
       MPH_ONFAILURE_GOTO9999(info)
       Call XMPH_dm_Create(bound_sol,bound_ndof,1,bound_ndof,info  ) 
       MPH_ONFAILURE_GOTO9999(info)
       Do i=1,bound_sol%m
          bound_sol%v(i) = XMPH_FLOATZERO
       End Do
    End If

    !---------------------------------------------------------
    ! [1.3] Allocate FGMRES workspace "work"
    !---------------------------------------------------------

    ! Set lwork value
    lwork = restrt**2 + restrt*(2*bound_ndof+5) + 5* bound_ndof + 1
     ! the workspace should be large enough to store the m dot-products 
    ! Allocate
    If ((fgmicntl(4).Eq.2).Or.(fgmicntl(4).Eq.3)) Then
       lwork = lwork + restrt
     Else
        lwork = lwork +1
     End If
    
    Allocate(work(lwork), STAT = info )
    CHCKALLOC(info)
    MPH_ONFAILURE_GOTO9999(info)
    
    !---------------------------------------------------------
    ! [1.4] Setup the data
    !---------------------------------------------------------
    ! Setup the rhs
    ! Work(1:bound_ndof)             = local initial guess
    ! Work(bound_ndof+1:2*bound_ndof) = local right-hand-side
    If (HasInitGuess) Call XMPH_ARITHcopy &
         (bound_ndof, bound_sol%v(1),1, work(1), 1)

    Call XMPH_ARITHcopy&
         (bound_ndof,bound_rhs%v(1),1,work(bound_ndof+1), 1)

    ! initialize vectors
    info = 0
    Call XMPH_dm_nullify(x,info)
    Call XMPH_dm_nullify(y,info)
    Call XMPH_dm_nullify(z,info)
    MPH_ONFAILURE_GOTO9999(info)

    x%ld = bound_ndof
    x%m  = bound_ndof
    x%n  = 1

    y%ld = bound_ndof
    y%m  = bound_ndof
    y%n  = 1

    z%ld = bound_ndof
    z%m  = bound_ndof
    z%n  = 1

    ! Associate aliases
    revcom => fgmirc(1)
    colx   => fgmirc(2)
    coly   => fgmirc(3)
    colz   => fgmirc(4)
    nbscal => fgmirc(5)

    !-----------------------------------------------------------------
    ! [2] Call the solver, Iterate until solution is found
    !-----------------------------------------------------------------
    fgmiinfo(2) = -1
    last_failed_it = -1
    fgmirc(8)=0
    revcom = 1 

    Do While ((revcom > 0).And.(info >= 0))

       !----------------------------------------------------
       ! [2.1] Call the driver
       !----------------------------------------------------


10       Call Drive_XMPH_ARITHfgmres(         &
            intrf_ndof,bound_ndof,  &
            restrt,lwork,work,    &
            fgmirc,fgmicntl,fgmrcntl, &
            fgmiinfo,fgmrinfo)

       !----------------------------------------------------
       ! [2.2] Associate arrays
       !----------------------------------------------------

       y%v => work(coly : (coly+bound_ndof-1) )
       z%v => work(colz : (colz+bound_ndof-1) )

       If ( revcom /= 4)Then
          x%n   = 1
          x%v => work(colx : (colx+bound_ndof-1) )
       Else
          x%n   =  nbscal
          x%v => work(colx : (colx + bound_ndof*nbscal -1) )
       End If



       !----------------------------------------------------
       ! [2.2.1] Fault simulation
       !----------------------------------------------------
       !check if the iteration will be faulty
       !and set flag to restart fgmres 
       If((fgmirc(8) .Eq. 0) .And. (fgmiinfo(2) .Gt. last_failed_it)) Then
          last_failed_it = fgmiinfo(2)
          Call XMPH_IsFault(mphs, fgmiinfo(2),fault_info)
          If(fault_info .Eq. FAULT_DOUBLE) Then !multiple fault
             fgmirc(8) = 1
          Else If(fault_info .Eq. FAULT_SIMPLE)Then !single fault
             Call XMPH_Check_Fault(mphs,lwork, work, fgmiinfo(2),fgmicntl)              
          End If
       End If

       ! Check if the fault can be taken into account, i.e., the current iterate has been                                                                                                              
       !computed and the internal variables have been reset to allow for a new call           

       If ((fgmirc(1).Eq.0).And.(fgmirc(8).Eq.1))Then !  .And.(fgmrinfo.Gt.fgmrcntl(1))) Then
          fgmicntl(5) = 1
          fgmirc(8) = 0
          Call XMPH_Check_Fault(mphs, lwork, work, fgmiinfo(2),fgmicntl) 
          Goto 10
       End If

       !----------------------------------------------------
       ! [2.3] Do what the driver asked to perform.
       !----------------------------------------------------
       !
       ! 1     - Perform Matrix Vector Product : z <-- Schur . x
       ! 2     - Apply left preconditioning    : z <-- Precond . x
       ! 3     - Apply right precondiotioning  : z <-- Precond . x
       ! 4     - Perform scalar product        : z <-- x * y
       ! 0     - Exit                          : Solving step ended
       ! other - Throw Error : unknown revcom
       !

       Select Case (revcom)
       Case (0); Continue                      ! Exit    
       Case (1); Call XMPH_SCHUR_SchurVectorProduct (mphs, x, z, info)
       Case (2); CHCKASSRT(.False., info )     ! MaPHyS does not support left preconditioning
       Case (3); Call XMPH_SCHUR_PcdVectorProduct   (mphs, x,    z, info)
       Case (4); Call XMPH_SCHUR_VectorScalarProduct(mphs, x, y, z, info)
       Case Default; CHCKASSRT(.False., info ) ! Undocumented revcom
       End Select

    End Do
            
    !---------------------------------------------------------------------------
    ! [-] Checking the outputs 
    !---------------------------------------------------------------------------    

    !-----------------------------------------------
    ! [--] Checking revcom, print a message on error
    !-----------------------------------------------

    Select Case (revcom)
    Case (0); Continue ! Gmres exited correctly
    Case (1); Call MPH_Log(MSG_ERROR,"Error occured in matrix/vector Product")
    Case (2); Call MPH_Log(MSG_ERROR,"MaPHyS does not support left preconditioning")
    Case (3); Call MPH_Log(MSG_ERROR,"Error occured in right preconditioning")
    Case (4); Call MPH_Log(MSG_ERROR,"Error occured in scalar product")
    Case Default; Call MPH_LogWithInfo(MSG_ERROR,revcom,&
            "From PackFGMRES: undocumented revcom code. revcom =")
    End Select
    MPH_ONFAILURE_GOTO9999(info)
    
    !----------------------------------------------------
    ! [--] Checking fgmiinfo(1), the return status of PackGMRES
    !----------------------------------------------------

    Select Case (fgmiinfo(1))
    Case ( 0) ! convergence achieved
       info = MPH_SUCCESS
       Continue
    Case (-4) ! non convergence. Consider it as a warning and not an error
       info = MPH_SUCCESS + 1
       Call MPH_LogWithInfo(MSG_WARNING,fgmicntl(7),&
            "From PackFGMRES: convergence not achieved. Nb iterations:")
    Case Default
       Call MPH_LogWithInfo(MSG_ERROR,fgmiinfo(1),"PackFGRMES exited with error code =")
       CHCKASSRT(.False., info )
    End Select
    MPH_ONFAILURE_GOTO9999(info)

    !-----------------------------------------------------------------
    ! [-] Exit routine
    !-----------------------------------------------------------------

    ! Save data
    Call XMPH_ARITHcopy(bound_ndof, work(1),1, bound_sol%v(1),1 )
    if (rank .eq. 0 ) Then    
       close(779) 
    End If
    mphs%IINFOG( IINFOG_ITS_NbIters ) = fgmiinfo(2)
    mphs%RINFOG( RINFOG_ITS_BckErr  ) = fgmrinfo

    mem = INT(lwork*XMPH_FLOATBYTESIZE,8)
    Call MPH_mem_add2usage(mphs%mem, mem )
    mphs%IINFO ( IINFO_ITS_MEMPEAK  ) = byte2Mbyte(mem)
    mphs%IINFO ( IINFO_ITS_MEMUSED  ) = byte2Mbyte(mem)

    ! Free memory
9999 Continue
    If (Associated(work)) Deallocate( work )

  End Subroutine XMPH_SCHUR_FGMRES_Solve



  ! [+] routine :  XMPH_SCHUR_CG_solve -----------------------------------------
  !
  !> solve the interface linear system using CG
  !!
  !! Solve the Schur complement system (distributed) with CG (Conjugate Gradient).
  !! where the system is :
  !!           S. x = rhs
  !!
  !!-----
  !!
  !! @param[in,out] mphs  
  !!
  !!       the maphys instance which contains :
  !!       - the parameters of the Iterative solver
  !!       - the Schur "S"      (which have multiple representation)
  !!       - Its Preconditioner (which have multiple representation)
  !!       - necessary data to communicate between processors.
  !!
  !! @param[in    ] bound_rhs  
  !!
  !!       it holds the local part of the right-hand-side "rhs"
  !!
  !! @param[in,out] bound_sol
  !!
  !!       - On input,  it may holds the local part of the initial guess x0
  !!       - On output, it holds the local part of the solution x
  !!
  !! @param[   out] info 
  !!       
  !!       Integer describing the routine status
  !!
  !!-----
  !!
  !! @note
  !! See also the documentation of PackCG at :
  !!   http://www.cerfacs.fr/algor/reports/2003/TR_PA_03_03.pdf
  !!
  !! @author Azzam Haidar
  !! @author Luc   Giraud
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_SCHUR_CG_solve(mphs, bound_rhs, bound_sol, info)

    !* Module used *!

    Use XMPH_maphys_type
    Use MPH_maphys_enum
    Use XMPH_schur_aux_mod
    Use XMPH_dense_matrix_mod
    Use mph_log_mod
    Use MPH_mem_mod
    Implicit None

    !* Arguments *!

    Type(XMPH_maphys_t)        , Intent(inout) :: mphs
    Type(XMPH_dense_matrix_t)  , Intent(in   ) :: bound_rhs 
    Type(XMPH_dense_matrix_t)  , Intent(inout) :: bound_sol 
    Integer               , Intent(  out) :: info

    !* Local Variables *!

    ! Scalars
    Integer    :: rank,master
    Integer    :: verb
    MPH_INT :: lwork
    MPH_INT :: intrf_ndof
    MPH_INT :: bound_ndof
    MPH_INT :: interior_ndof
    MPH_INT :: domain_ndof
    Integer    :: i
    Logical    :: HasInitGuess
    Integer(kind=8) :: mem

    ! Arrays
    Integer, Target :: cgirc   (5)
    Integer         :: cgicntl (7)
    Real(Kind=XMPH_FLOATKIND) :: cgrcntl (3)

    Integer         :: cgiinfo (3)
    Real(Kind=XMPH_FLOATKIND) :: cgrinfo (3)

    XMPH_FLOAT, Pointer :: work     (:)

    ! Aliases (Pointers)

    ! to cgirc (reverse communication)
    Integer, Pointer :: revcom
    Integer, Pointer :: colx
    Integer, Pointer :: coly
    Integer, Pointer :: colz

    ! to cg workspace "work"
    Type(XMPH_dense_matrix_t) :: x 
    Type(XMPH_dense_matrix_t) :: y 
    Type(XMPH_dense_matrix_t) :: z

    !- End of header -------------------------------------------------

    !-----------------------------------------------------------------
    ! [1] Initialize local variables / CG package
    !-----------------------------------------------------------------

    !---------------------------------------------------------
    ! [1.1] Get differents sizes & options
    !---------------------------------------------------------

    intrf_ndof    = mphs%lc_domain%gballintrf
    bound_ndof    = mphs%lc_domain%myndofintrf
    interior_ndof = mphs%lc_domain%myndofinterior
    domain_ndof   = mphs%lc_domain%myndof
    HasInitGuess = (mphs%IKEEP(IKEEP_ITS_InitGuess) == 1)

    !---------------------------------------------------------
    ! [1.2] Set CG options
    !---------------------------------------------------------


    ! Get default options
    Call init_XMPH_ARITHcg(cgicntl,cgrcntl)

    ! Set specific options
    cgrcntl(1) = Real( mphs%rkeep(RKEEP_ITS_Tolerance),KIND=XMPH_FLOATKIND)

    cgicntl(1) = 6 ! print errors 
    cgicntl(2) = 0 ! no warning messages
    cgicntl(3) = 0 ! no convergence history
    If ( mphs%ikeep(IKEEP_Pcd_Strategy ) == PCD_STRATEGY_isNone )&
         cgicntl(4) = 0 ! no preconditioner
    If ( mphs%ikeep(IKEEP_Pcd_Strategy ) /= PCD_STRATEGY_isNone )&
         cgicntl(4) = 1 ! left preconditioner
    cgicntl(5) = mphs%ikeep(IKEEP_ITS_InitGuess   ) ! initial guess
    cgicntl(6) = mphs%ikeep(IKEEP_ITS_MaxIters    ) ! maximum number of iterations
    cgicntl(7) = 0 ! no eigenvalues

    ! Set output 
    rank= mphs%ikeep(IKEEP_MPIRANK)
    master= mphs%ikeep(IKEEP_HOSTRANK)
    Call mph_log_GetVerbosity(verb)
    If ( (verb >= MSG_STD).and.(rank == master) )Then
       ! print warning messages
       Call mph_log_GetUnit(MSG_WARNING, cgicntl(2))
       ! print convergence history to stdout
       Call mph_log_GetUnit(MSG_STD    , cgicntl(3))
    End If

    ! Create solution if no init guess
    If (HasInitGuess .eqv. .False.) Then

       Call XMPH_dm_Nullify(bound_sol,info) 
       MPH_ONFAILURE_GOTO9999(info)

       Call XMPH_dm_Create(bound_sol,bound_ndof,1,bound_ndof,info  ) 
       MPH_ONFAILURE_GOTO9999(info)

       Do i=1,bound_sol%m
          bound_sol%v(i) = XMPH_FLOATZERO
       End Do

    End If

    !---------------------------------------------------------
    ! [1.3] Allocate CG workspace "work"
    !---------------------------------------------------------

    ! Set lwork value
    lwork = 6*intrf_ndof + 1
    If (cgicntl(7) /= 0) lwork = lwork + 2*(cgicntl(6)+1)

    ! Allocate
    Allocate(work(lwork), STAT = info )
    CHCKALLOC(info)
    MPH_ONFAILURE_GOTO9999(info)

    !---------------------------------------------------------
    ! [1.4] Setup the data
    !---------------------------------------------------------


    ! Setup the rhs
    ! Work(1:bound_ndof)             = local initial guess
    ! Work(bound_ndof+1:2*bound_ndof) = local right-hand-side
    If (HasInitGuess) Call XMPH_ARITHcopy &
         (bound_ndof, bound_sol%v(1),1, work(1), 1)

    Call XMPH_ARITHcopy &
         (bound_ndof,bound_rhs%v(1),1,work(bound_ndof+1),1)

    ! initialize vectors
    Call XMPH_dm_nullify(x,info)
    MPH_ONFAILURE_GOTO9999(info)

    Call XMPH_dm_nullify(y,info)
    MPH_ONFAILURE_GOTO9999(info)

    Call XMPH_dm_nullify(z,info)
    MPH_ONFAILURE_GOTO9999(info)

    x%ld = bound_ndof
    x%m  = bound_ndof
    x%n  = 1

    y%ld = bound_ndof
    y%m  = bound_ndof
    y%n  = 1

    z%ld = bound_ndof
    z%m  = bound_ndof
    z%n  = 1

    ! Associate aliases
    revcom => cgirc(1)
    colx   => cgirc(2)
    coly   => cgirc(3)
    colz   => cgirc(4)

    !-----------------------------------------------------------------
    ! [2] Call the solver, Iterate until solution is found
    !-----------------------------------------------------------------

    revcom = 1 
    Do While ((revcom > 0).And.(info >= 0))

       !----------------------------------------------------
       ! [2.1] Call the driver
       !----------------------------------------------------

       Call drive_XMPH_ARITHcg(         &
            intrf_ndof,bound_ndof,  &
            lwork,work,    &
            cgirc,cgicntl,cgrcntl, &
            cgiinfo,cgrinfo)

       !----------------------------------------------------
       ! [2.2] Associate arrays
       !----------------------------------------------------

       x%v => work(colx : (colx+bound_ndof-1) )
       y%v => work(coly : (coly+bound_ndof-1) )
       z%v => work(colz : (colz+bound_ndof-1) )

       !----------------------------------------------------
       ! [2.3] Do what the driver asked to perform.
       !----------------------------------------------------
       !
       ! 1     - Perform Matrix Vector Product : z <-- Schur . x
       ! 2     - Apply right precondiotioning  : z <-- Precond . x
       ! 3     - Perform scalar product        : z <-- x * y
       ! 0     - Exit                          : Solving step ended
       ! other - Throw Error : unknown revcom
       !

       Select Case (revcom)
       Case (0); Continue
       Case (1); Call XMPH_SCHUR_SchurVectorProduct (mphs, x,    z, info)
       Case (2); Call XMPH_SCHUR_PcdVectorProduct   (mphs, x,    z, info)
       Case (3); Call XMPH_SCHUR_VectorScalarProduct(mphs, x, y, z, info)
       Case DEFAULT ; CHCKASSRT(.False.,info)
       End Select

    End Do

    !---------------------------------------------------------------------------
    ! [-] Checking the outputs 
    !---------------------------------------------------------------------------    

    !-----------------------------------------------
    ! [--] Checking revcom, print a message on error
    !-----------------------------------------------

    Select Case (revcom)
    Case (0); Continue ! Gmres exited correctly
    Case (1); Call MPH_Log(MSG_ERROR,"Error occured in matrix/vector Product")
    Case (2); Call MPH_Log(MSG_ERROR,"Error occured in preconditioning")
    Case (3); Call MPH_Log(MSG_ERROR,"Error occured in scalar product")
    Case Default; Call MPH_LogWithInfo(MSG_ERROR,revcom,&
            "From PackCG: undocumented revcom code. revcom =")
    End Select
    MPH_ONFAILURE_GOTO9999(info)
    
    !----------------------------------------------------
    ! [--] Checking gmiinfo(1), the return status of PackGMRES
    !----------------------------------------------------

    Select Case (cgiinfo(1))
    Case ( 0) ! convergence achieved
       info = MPH_SUCCESS
       Continue
    Case (-3) ! non convergence. Consider it as a warning and not an error
       info = MPH_SUCCESS + 1
       Call MPH_LogWithInfo(MSG_WARNING,cgicntl(6),&
            "From PackCG: convergence not achieved. Nb iterations:")
    Case Default
       Call MPH_LogWithInfo(MSG_ERROR,cgiinfo(1),"PackCG exited with error code =")
       CHCKASSRT(.False., info )
    End Select
    MPH_ONFAILURE_GOTO9999(info)

    !-----------------------------------------------------------------
    ! [3] Exit routine
    !-----------------------------------------------------------------

    ! Save data
    Call XMPH_ARITHcopy(bound_ndof, work(1),1, bound_sol%v(1),1 )
    mphs%IINFOG( IINFOG_ITS_NbIters ) = cgiinfo(2)
    mphs%RINFOG( RINFOG_ITS_BckErr  ) = cgrinfo(1)

    mem = INT(lwork*XMPH_FLOATBYTESIZE,8)
    Call MPH_mem_add2usage(mphs%mem, mem )
    mphs%IINFO ( IINFO_ITS_MEMPEAK  ) = byte2Mbyte(mem)
    mphs%IINFO ( IINFO_ITS_MEMUSED  ) = byte2Mbyte(mem)

    ! Free memory
9999 Continue
    If (Associated(work)) Deallocate( work )

  End Subroutine XMPH_SCHUR_CG_Solve



End Module XMPH_SCHUR_mod


