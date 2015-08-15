! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"
!! module for fault pluging.
!!
!! 
Module XMPH_fault_mod

  !* Modules * !
  Use mph_error_mod
  Implicit None
  !* Private constants *!
  Character(len=MAPHYS_STRL), Private, Parameter :: &
       FLNAME = "xmph_fault_mod.F90"

  !* List of routines *!
  Public   :: XMPH_CHECK_FAULT
  Public   :: XMPH_ISFAULT
Contains
  
  Subroutine XMPH_Check_Fault( &! intents
       mphs,                   & ! in
       lwork,                  & ! inout
       work ,                  & ! in
       iter,                   & ! in
       gmicntl                 & ! in
       )
    
    !*modules *!
    Use XMPH_maphys_type
    Use MPH_maphys_enum
    Use XMPH_dense_matrix_mod  
    Use MPH_time_mod
!!DBG
    Use XMPH_dds_mod
    Use XMPH_schur_aux_mod
    Implicit None 
    Include 'mpif.h'
    !* Arguments *!
    Type(XMPH_maphys_t)   , Intent(inout) :: mphs
    Integer               , Intent(in )   :: lwork
    XMPH_FLOAT            :: work(*)
    Integer               , Intent(in )   :: iter
    Integer               , Intent(in )  :: gmicntl(*)
    !* Local variables *!
    Integer      :: i, FT_strategy 
    Integer      :: Diff
    Integer      :: myrank
    Integer      :: bound_ndof
    Integer      :: restrt ! restart parameter
    Integer      :: info
    Integer      :: ArrayOfFaultyPorc(3)
    Integer      :: FailedNeighbs(3)
    Integer      :: intrf_ndof
    Integer      :: NbOfFailedNeighb
    Integer      :: NbOfFaultyProc
    Integer      :: Iamfailed
    XMPH_FLOAT   :: test

    Type(XMPH_dense_matrix_t)  :: bound_rhs 
    Type(XMPH_dense_matrix_t)  :: x
    !* initialize variables*!
    Diff = 0
    Iamfailed = 0
    FT_strategy = mphs%icntl(29)
    myrank= mphs%ikeep(IKEEP_MPIRANK)
    bound_ndof    = mphs%lc_domain%myndofintrf
    restrt  = mphs%ikeep(IKEEP_ITS_Restart) 
    intrf_ndof    = mphs%lc_domain%gballintrf        
    If (restrt .gt. intrf_ndof) restrt = intrf_ndof
    
    !---------------------------------------------------------
    ! [1] Fault detection
    !---------------------------------------------------------
    ! [1.1] Timers
    !sett initialize timer
    Call MPH_time_start(mphs%rinfo (RINFO_TIMING_Fault))
    Call XMPH_GenerateFault(mphs,ArrayOfFaultyPorc,NbOfFaultyProc,Iamfailed,iter)
    Call MPH_time_stop(mphs%rinfo (RINFO_TIMING_Fault))
    Call MPH_time_start(mphs%rinfo (RINFO_TIMING_Interp))
    !---------------------------------------------------------
    ! [2] Fault recovery
    !---------------------------------------------------------
    If (NbOfFaultyProc .Eq. 1)  Then 
       !---------------------------------------------------------
       ! [2.1] single fault case
       !---------------------------------------------------------
       If (Iamfailed .Eq. 1) Then 
          !---------------------------------------------------------
          ! [2.1.a] set  workspace to zero
          !---------------------------------------------------------
          Do i =1,lwork
             work(i) = XMPH_FLOATZERO
          End Do
          
          !---------------------------------------------------------
          ! [2.1.b] Receive data from neighbors
          !---------------------------------------------------------   
          Call XMPH_FRecv_Neighbors(mphs, lwork, work, NbOfFaultyProc, gmicntl, iter, info)
       Else
          !---------------------------------------------------------
          ! [2.1.c] Send data to failed neighbor
          !---------------------------------------------------------           
          Call XMPH_IsNeighbor(mphs, ArrayOfFaultyPorc, NbOfFaultyProc,&
               FailedNeighbs, NbOfFailedNeighb,info)
          If(info == 1 ) Then  ! The faulty processor is my neighbor
             Call XMPH_FSend_Neighbors(mphs, lwork, work,&
                  FailedNeighbs,NbOfFailedNeighb, gmicntl, iter,info)           
          End If
       End If
    End if
    
    !---------------------------------------------------------
    ! [2.2] multiple faults case
    !---------------------------------------------------------
    
    If(NbOfFaultyProc .gt. 1) Then  

       If (Iamfailed .Eq. 1) Then
          !-------------------------------------------------
          ! [2.2.a] Set  workspace to zero except the rhs 
          !---------------------------------------------------------
          If(FT_strategy .Le. 1) Then
             Do i =1,bound_ndof
                work(i) = XMPH_FLOATZERO
             End Do
          End If
          Do i=2*bound_ndof+1,lwork
            work(i) = XMPH_FLOATZERO
          End DO
       End If
       
       Call XMPH_IsNeighbor(mphs, ArrayOfFaultyPorc, NbOfFaultyProc,&
            FailedNeighbs, NbOfFailedNeighb,info)
       
       ! Restore lost date and rhs contribution
       If ((info .Eq. 1) .Or. (Iamfailed .Eq. 1)) Then
           If( (FT_strategy .EQ. 0) .Or. (FT_strategy .EQ. 2) ) Then
              Call XMPH_RestoreX(mphs, work,FailedNeighbs,NbOfFailedNeighb,Iamfailed)
           Elseif(FT_strategy .EQ. 1) Then 
              Call XMPH_ExchContrib(mphs, work,FailedNeighbs,NbOfFailedNeighb,Iamfailed)
           End If
       End If
       ! INTERPOLATION
       If((Iamfailed .Eq. 1) .AND. (NbOfFailedNeighb .Gt. 0) &
            .AND. (FT_strategy .Eq. 1)) Then
          Call XMPH_Interpolate(mphs,work, ArrayOfFaultyPorc, NbOfFaultyProc,&
               FailedNeighbs, NbOfFailedNeighb,info)
       End If
    End If
    ![--] timers                                                                                                                   
    Call MPH_time_stop(mphs%rinfo (RINFO_TIMING_Interp))
    
    !---------------------------------------------------------
    ! [4] Check the correctness of restored data
    !---------------------------------------------------------
    If((NbOfFaultyProc .gt. 0) .And. (Iamfailed .Eq. 1)) Then
       !***DBG

       Diff = 0
    End If
    info =0
    !----------------------------------------------------------------------
    ! [5] Exit routine
    !----------------------------------------------------------------------
    
9999 Continue
       
  END SUBROUTINE XMPH_Check_Fault
  
     !XMPH_ISFAULT return
     ! 0 if there is no fault
     !else the number of faults
     SUBROUTINE XMPH_ISFAULT(mphs, iter,info)
       
       Use XMPH_maphys_type
       Use MPH_maphys_enum
       Implicit None
       !*Arguments*!
       Type(XMPH_maphys_t) ,Intent(in ) :: mphs
       Integer             ,Intent(in ) :: iter
       Integer             ,Intent(inout) :: info

       Integer FT_taux
       Integer FT_MaxFault
       Integer nbfault
       FT_taux = mphs%icntl(28)
       FT_MaxFault = mphs%icntl(33)
       nbfault = mphs%icntl(34)
       info = 0
       If((mod(iter, FT_taux) .eq. 0) .AND. ((FT_MaxFault*FT_taux) .ge. iter))  info = nbfault
       
     END SUBROUTINE XMPH_ISFAULT
     

     SUBROUTINE XMPH_GenerateFault(mphs,ArrayOfFaultyPorc, NbOfFaultyProc,Iamfailed,iter)

       
       Use XMPH_maphys_type
       Use MPH_maphys_enum
       
       Implicit None
       Include'mpif.h'
       
       !Arguments *!
       Type(XMPH_maphys_t) ,Intent(in ) :: mphs
       Integer             ,Intent(inout) :: ArrayOfFaultyPorc(3)
       Integer             ,Intent(out )  :: nbOfFaultyProc
       Integer             ,Intent(in  )  :: iter
       Integer             ,Intent(inout ):: Iamfailed
       !local variables*!
       Integer i, j, nNeighb
       Integer myrank, root
       Integer nbproc
       Integer minCriterion
       Integer candidate
       Integer info
       Integer Faulty(2)
       Integer nbfault
       Integer k
       Integer Min_k
       Integer Best_candidate
       Integer bound_ndof
       Real    rand_number, aux(iter)
       Integer, Pointer      :: failedproc(:)
       Integer, Pointer      :: failedcandidate(:)
       Integer, Pointer :: indexVi       (:)
       Integer, Pointer :: ptr_Index_Intrfc (:)
       
       nbproc= mphs%ikeep(IKEEP_NBDOMAINS)
       bound_ndof    = mphs%lc_domain%myndofintrf
       myrank= mphs%ikeep(IKEEP_MPIRANK)
       nbfault = mphs%icntl(34)
       nNeighb = mphs%lc_domain%mynbvi
       indexVi => mphs%lc_domain%myindexVi
       ptr_Index_Intrfc => mphs%lc_domain%myptrindexVi 
       root = 0
       candidate = 0
       Min_k=10000000
       !* memory allocation
       Allocate(failedproc(nbproc), STAT = info)
       CHCKASSRT( info == 0, info )
       If( info < 0 ) Goto 9999
       
       !* memory allocation
       Allocate(failedcandidate(nbproc), STAT = info)
       CHCKASSRT( info == 0, info )
       If( info < 0 ) Goto 9999

       
       NbOfFaultyProc = 0
       !* All processor determine the one with less neighbors 
       Call MPI_Allreduce(nNeighb, minCriterion, 1, MPI_INT, MPI_MIN, mphs%comm, info)
       
       If (minCriterion .Eq. nNeighb) Then
          candidate = 1
       End If
       
       Call MPI_Gather(candidate, 1, MPI_INT, &
            failedcandidate, 1, MPI_INT, root, mphs%comm, info); 
       ASSRT( info == MPI_SUCCESS )

       !* The proc of rank 0 choose the first faulty 
       !* Proc ramdomly
       If (myrank .Eq. 0) Then
          ! Select the processor with smallest rank
          Do i = 1, nbproc
             If (failedcandidate(i) .Eq. 1) Then
                Faulty(1) = i-1
                !Faulty(1) = nbproc/2
             End If
          End Do
       End If
       !* Broadcast the id of the faulty proc 
       Call MPI_Bcast(Faulty(1), 1, MPI_INT, root, mphs%comm, info)
       ASSRT( info == MPI_SUCCESS )
       !*Each proc If he is the failed proc
       If( Faulty(1) .eq. myrank) Then
          Iamfailed = 1
       End If

       ! determine the second failed proc if 
       ! we simulate multiple faults
       If( nbfault .Eq. FAULT_DOUBLE) Then
          !Send all number of neighbor to failed one
          Call MPI_Gather(nNeighb, 1, MPI_INT, &
               failedcandidate, 1, MPI_INT, Faulty(1), mphs%comm, info); 
          ASSRT( info == MPI_SUCCESS )
          If (Iamfailed .Eq. 1) Then
             If (nNeighb .Gt. 0) Then
                Best_candidate =1
                Do i= 1, nNeighb
                   k=  failedcandidate(indexVi(i)+1)
                   If (k .le. Min_k ) Then
                      Min_k = k
                      Best_candidate = i
                   End If
                End Do
                Faulty(2) =  indexVi(Best_candidate)
             End If
          End If
          
          !* Broadcast the id of the second faulty proc
          root = Faulty(1)
          Call MPI_Bcast(Faulty(2), 1, MPI_INT, root, mphs%comm, info)          
          ASSRT( info == MPI_SUCCESS )

          !*Each proc check If it is the failed proc
          If( Faulty(2) .eq. myrank) Then
             Iamfailed = 1
          End If
       End If
       ! Each processor send Iamfailed
       ! to others 
       Call MPI_Allgather(Iamfailed, 1, MPI_INT, &
            failedproc, 1, MPI_INT, mphs%comm, info); 
       ASSRT( info == MPI_SUCCESS )
       
       
       ! count number of fault
       Do i = 1, nbproc
          If (failedproc(i) .Eq. 1) Then
             NbOfFaultyProc = NbOfFaultyProc+1
             ArrayOfFaultyPorc(NbOfFaultyProc) = i-1
          End If
       End Do
       
9999 continue
    !Free memory
    If (Associated(failedproc)) Deallocate(failedproc)

  End SUBROUTINE XMPH_GENERATEFAULT


  !* MPH_IsNeighbor check if the faulty proc
  !* is my neighbor
  SUBROUTINE XMPH_IsNeighbor( & ! intents
       mphs,                 & ! in 
       ArrayOfFaultyPorc,    & ! in
       NbOfFaultyProc,       & ! in
       FailedNeighbs,        & ! in
       NbOfFailedNeighb,     & ! in
       info                  & ! out 
       )
       
    !* Module(s) *!
    Use XMPH_maphys_type
    Use mph_error_mod        
    Use MPH_maphys_enum
    
    Implicit None
    Include 'mpif.h'
    
    !* Arguments *!
    Type(XMPH_maphys_t)   , Intent(in) :: mphs
    Integer               , Intent(in  )  :: ArrayOfFaultyPorc(3)
    Integer               , Intent(in  )  :: NbOfFaultyProc
    Integer               , Intent(out  )  :: FailedNeighbs(3)
    Integer               , Intent(out )  :: NbOfFailedNeighb
    Integer               , Intent(out ) :: info    
    !* Local variables *!
    
    ! Scalars
    Integer :: nNeighb 
    Integer :: neighb  
    Integer :: failedrank
    Integer :: i,myrank
    ! Arrays 
    Integer     , Pointer :: indexVi       (:)
    
    !initialization
    info = 0
    myrank= mphs%ikeep(IKEEP_MPIRANK)
    nNeighb      = mphs%lc_domain%mynbvi
    indexVi      => mphs%lc_domain%myindexVi
    NbOfFailedNeighb = 0
    Do i =1,NbOfFaultyProc
       failedrank = ArrayOfFaultyPorc(i) 
       Do neighb=1, nNeighb
          If ((failedrank .Eq. indexVi(neighb)) .And. (failedrank .Ne. myrank)) Then
             info = 1
             NbOfFailedNeighb = NbOfFailedNeighb +1
             FailedNeighbs(NbOfFailedNeighb) = failedrank
          End If
       End Do
    End DO
9999 Continue
      
  END SUBROUTINE XMPH_IsNeighbor
  
  
  SUBROUTINE   XMPH_Send_Neighbors( & ! intents
       mphs,                 & ! inout
       work,                 & ! in
       gmicntl,              & ! in 
       failediter,           & ! in 
       info                  & ! out
       )       
       
    !* Module(s) *!
    Use XMPH_maphys_type
    Use mph_error_mod        
    Use MPH_maphys_enum
    Implicit None
    Include 'mpif.h'
    
    !* Arguments *!
    Type(XMPH_maphys_t)    , Intent(inout) :: mphs
    XMPH_FLOAT                  :: work(*)
    Integer                , Intent( in ) :: gmicntl(8)
    Integer                , Intent(in )  :: failediter
    Integer                , Intent(out )  :: info
    
    !* Local variables *!

    ! Scalars
    Integer, Parameter :: PresetMPITagg = 77
    Integer :: tagg, myrank 
    Integer :: nNeighb 
    Integer :: neighb  
    Integer :: indexFailed_Neighbor
    Integer :: mpiReq
    Integer :: restrt ! restart parameter
    Integer :: bound_ndof  ! ndof of the interface
    Integer :: intrf_ndof
    Integer :: sizedot
    Integer :: failedrank
    MPH_INT :: maxSizeIntrf
    MPH_INT ::  sendBuffSize
    MPH_INT ::  sendBuffStart
    MPH_INT ::  sendBuffend
    MPH_INT :: i,j, ij ,k 
    ! Arrays 
    Integer     , Pointer :: indexVi       (:)
    Integer     , Pointer :: ptr_Index_Intrfc (:)
    Integer     , Pointer :: Index_Intrfc(:)
    XMPH_FLOAT  , Pointer :: sendbuff             (:)
    Integer :: mpistatus(MPI_STATUS_SIZE) 
    
    !-------------------------------------------------------------------
    ![1] Initialize local arguments
    !-------------------------------------------------------------------
    ! restart parameter
    restrt  = mphs%ikeep(IKEEP_ITS_Restart) 
    ! communication  
    info = 0
    myrank= mphs%ikeep(IKEEP_MPIRANK)
    tagg  = PresetMPITagg
    sendBuffSize = 0    
    ! Domain related
    nNeighb          = mphs%lc_domain%mynbvi
    indexVi          => mphs%lc_domain%myindexVi
    ptr_Index_Intrfc => mphs%lc_domain%myptrindexVi 
    indexFailed_Neighbor = -1
    bound_ndof    = mphs%lc_domain%myndofintrf
    Index_Intrfc     =>mphs%lc_domain%myindexIntrf
    intrf_ndof    = mphs%lc_domain%gballintrf    
    failedrank = 2
    !! check restrt
    If (restrt .gt. intrf_ndof) restrt = intrf_ndof
    If ((gmicntl(5).Eq.2).Or.(gmicntl(5).Eq.3)) Then
       sizedot = restrt
    Else 
       sizedot = 1
    End If

    !-------------------------------------------------------------------------
    ! [2] Determination of the interface with the failed neighbor 
    !-------------------------------------------------------------------------
    
    ! search the local index of
    ! the failed neighbor 
    Do neighb=1, nNeighb
       If (failedrank .Eq. indexVi(neighb)) Then
          indexFailed_Neighbor = neighb
          Exit  
       End If
    End Do

    ! k is the size of the interface, 
    k = ptr_Index_Intrfc(indexFailed_Neighbor+1) &
         - ptr_Index_Intrfc(indexFailed_Neighbor)
    
    ! compute size of the send buffer
    If(gmicntl(8) .Eq. 1) Then 
       sendBuffSize = (restrt+1)*(failediter+1) + (failediter+1)*k&
            + 5*k + 3*restrt + sizedot + 1
    Else
       sendBuffSize =(restrt+1)*(failediter+1) + (failediter+1)*k&
            + 6*k + 3*restrt + sizedot + 1
    End If
    ! Perform allocation
    Allocate( sendbuff(sendBuffSize), STAT = info)
    CHCKASSRT( info == 0, info )
    If( info < 0 ) Goto 9999



    !-------------------------------------------------------------------------
    ! [3] Copy and send data to the failed processor 
    !-------------------------------------------------------------------------
    
    !-------------------------------------------------------------------------
    ! [3.1] Copy distributed data (x,b, r, w,V) 
    !------------------------------------------------------------------------
  
    k = 1
    DO i= 1,(4+(failediter+1))
       Do j=ptr_Index_Intrfc(indexFailed_Neighbor),&
            ptr_Index_Intrfc(indexFailed_Neighbor+1)-1
          ij = index_intrfc(j) + bound_ndof*(i-1)
          sendbuff(k) = work(ij)
          k = k+1
       End Do
    End Do

    !-------------------------------------------------------------------------   
    ! [3.2] Copy local shared data (H,dot,ycurrent)
    !-------------------------------------------------------------------------
    If(gmicntl(8) .Eq. 1) Then
       sendBuffStart  = (4+restrt)*bound_ndof+1 !
    Else
       sendBuffStart = (5+restrt)*bound_ndof +1
    End If

    !copy failediter+1 vectors
    Do i=sendBuffStart, sendBuffStart+& 
         (restrt+1)*(failediter+1)-1
       sendbuff(k) = work(i)
       k= k+1
    End Do
    !copy last vector
    Do i= sendBuffStart +(restrt+1)*(restrt),&
         sendBuffStart +(restrt+1)*(restrt+1)-1
       sendbuff(k) = work(i)
       k= k+1
    End DO
    
    ! Copy dot
    sendBuffStart = sendBuffStart + (restrt+1)*(restrt+1)
    Do i=1,sizedot
       sendbuff(k)= work(sendBuffStart+i-1)
       k = k +1
    End Do
    !copy ycurrent
    sendBuffStart = sendBuffStart+ sizedot
    Do i=1,failediter+1
       sendbuff(k) = work(sendBuffStart + i - 1)
       k = k+1
    End DO
    !-------------------------------------------------------------------------
    ! [3.3] copy distributed data Xcurrent
    !-------------------------------------------------------------------------

    sendBuffStart = sendBuffStart + restrt
    sendBuffend = sendBuffStart + bound_ndof-1
    Do j=ptr_Index_Intrfc(indexFailed_Neighbor),&
         ptr_Index_Intrfc(indexFailed_Neighbor+1)-1
       
       sendbuff(k) = work(index_intrfc(j)+sendBuffStart-1)
       k = k+1
    End Do
    !-------------------------------------------------------------------------    
    ! [3.4] copy local shared data rotsin, rotcos
    !-------------------------------------------------------------------------
  
    sendBuffStart = sendBuffend 
    Do i = 1,2 
       DO j =1,failediter+1
          ij = j + sendBuffStart+ (i-1)*restrt
          sendbuff(k) = work(ij)
          k = k+1
       End Do
    End Do
    !-------------------------------------------------------------------------
    ! [4 ] send data
    !-------------------------------------------------------------------------
    CALL MPI_send(                    &
         sendbuff(1),sendBuffSize,XMPH_FLOATMPI,failedrank,&
         tagg,mphs%comm,info)
    ASSRT( info == MPI_SUCCESS )
    
  info = 0
  !----------------------------------------------------------------------
  ! [5] Exit routine
  !----------------------------------------------------------------------
  
9999 Continue
  !Free memory
  If (Associated(sendbuff)) Deallocate (sendbuff)

End SUBROUTINE XMPH_Send_Neighbors



SUBROUTINE XMPH_Recv_Neighbors(& ! intents
     mphs,                 & ! in
     work,                 & ! inout
     gmicntl,              & ! in
     failediter,           & ! in
     info                  & ! out
     )
     
  !* Module(s) *!
  Use XMPH_maphys_type
  Use mph_error_mod        
  Use MPH_maphys_enum
  Implicit None
  Include 'mpif.h'
  
  !* Arguments *!
  Type(XMPH_maphys_t)   , Intent(inout) :: mphs
  XMPH_FLOAT            , Intent(inout) :: work(*)
  Integer               , Intent(in )   :: gmicntl(8)
  Integer               , Intent(in )   :: failediter
  Integer               , Intent(out )  :: info
  
  !* Local variables *!
  ! Scalars
  Integer, Parameter :: PresetMPITagg = 77
  Integer :: tagg, myrank 
  Integer :: nNeighb 
  Integer :: neighb  
  Integer :: restrt ! restart parameter
  Integer :: bound_ndof  ! ndof of the interface
  Integer :: intrf_ndof
  Integer :: first
  Integer :: sizedot
  MPH_INT :: maxBuffSize
  MPH_INT :: maxInterSize
  MPH_INT ::  recvBuffSize
  MPH_INT ::  recvBuffStart
  MPH_INT ::  recvBuffend
  MPH_INT :: i,j, ij ,k,l 
  ! Arrays 
  Integer     , Pointer :: indexVi       (:)
  Integer     , Pointer :: ptr_Index_Intrfc (:)
  Integer     , Pointer :: Index_Intrfc     (:)
  XMPH_FLOAT  , Pointer :: recvbuff          (:)
  Integer     , Pointer :: mpiReq (:)
  Integer :: mpistatus(MPI_STATUS_SIZE) 
  
  !-------------------------------------------------------------------
  ![1] Initialize local arguments
  !-------------------------------------------------------------------
  ! restart parameter
  restrt  = mphs%ikeep(IKEEP_ITS_Restart) 
  ! communication  
  info = 0
  myrank= mphs%ikeep(IKEEP_MPIRANK)
  tagg  = PresetMPITagg
  recvBuffSize = 0    
  ! Domain related
  nNeighb          = mphs%lc_domain%mynbvi
  indexVi          => mphs%lc_domain%myindexVi
  ptr_Index_Intrfc => mphs%lc_domain%myptrindexVi 
  bound_ndof    = mphs%lc_domain%myndofintrf
  Index_Intrfc     => mphs%lc_domain%myindexIntrf
  intrf_ndof    = mphs%lc_domain%gballintrf
  maxInterSize = 0
  first = 1
   Allocate ( mpiReq(nNeighb), STAT = info )
    CHCKASSRT( info == 0, info )
     If( info < 0 ) Goto 9999
    !! check restrt
    If (restrt .ge. intrf_ndof)restrt = intrf_ndof
    If ((gmicntl(5).Eq.2).Or.(gmicntl(5).Eq.3)) Then
       sizedot = restrt
    Else 
       sizedot = 1
    End If
  !-------------------------------------------------------------------------
  ! [2] Detect message from neighbor and allocate receive buffer 
  !-------------------------------------------------------------------------
  
     !search maximum of interfaces size
     Do neighb=1,nNeighb 
        k = ptr_Index_Intrfc(neighb+1)-ptr_Index_Intrfc(neighb)
        If (k > maxInterSize) Then
           maxInterSize = k
        End If
     End Do
     
     If ( gmicntl(8) .Eq. 1) Then  
        maxBuffSize = restrt**2 + restrt*(maxInterSize+5)&
             + 5*maxInterSize +sizedot +  1
     Else
        maxBuffSize = restrt**2 + restrt*(maxInterSize+5)&
             + 6*maxInterSize + sizedot + 1
     End If
     ! Perform allocation
        Allocate (recvbuff(maxBuffSize),STAT = info )
        CHCKASSRT( info == 0, info )
        If( info < 0 ) Goto 9999

  Do l=1,nNeighb
     ! Get a message
     Call MPI_Probe(MPI_ANY_SOURCE,tagg,mphs%comm,mpiStatus,info)
     ASSRT( info == MPI_SUCCESS )
     
     Call MPI_Get_count(mpiStatus,XMPH_FLOATMPI,recvBuffsize,info)
     ASSRT( info == MPI_SUCCESS )

     ! Detect which neighbor this message corresponds
     Do neighb=1,nNeighb
	If ( indexVi(neighb)  ==  mpiStatus(MPI_SOURCE) ) Exit 
     End Do

     !-------------------------------------------------------------------------
     ! [3] receive and restore data 
     !-------------------------------------------------------------------------
     Call MPI_recv(recvbuff(1),recvBuffSize,&
          XMPH_FLOATMPI,mpiStatus(MPI_SOURCE),&
          tagg,mphs%comm,mpiReq(l),info)
     ASSRT( info == MPI_SUCCESS )
  
     !-------------------------------------------------------------------------
     ! [3.1] restore distributed data (x,b, r, w,V) 
     !-------------------------------------------------------------------------
     
     k = 1         
     Do i=1,(4+(failediter+1)) 
        Do j=ptr_Index_Intrfc(neighb),&
             ptr_Index_Intrfc(neighb+1)-1
           ij = index_intrfc(j) + bound_ndof*(i-1)
           work(ij) = recvbuff(k)
           k = k +1
        End Do
     End Do

     !-------------------------------------------------------------------------   
     ! [3.2] restore local shared data (H,dot,ycurrent)  
     !-------------------------------------------------------------------------
     
     If (gmicntl(8) .Eq. 1)Then
        recvBuffStart  = (4+restrt)*bound_ndof+1 !
     Else
        recvBuffStart = (5+restrt)*bound_ndof+1
     End If

     ! restore (failediter+1) vectors 
     Do i=recvBuffStart, recvBuffStart+&
          (restrt+1)*(failediter+1)-1
        if(first == 1) Then   
           work(i) =   recvbuff(k)  
        End If
        k = k+1
     End Do
     ! restore last vector of H
     Do i= recvBuffStart + (restrt+1)*(restrt),&
          recvBuffStart +(restrt+1)*(restrt+1)-1
        if(first == 1) Then   
           work(i) =   recvbuff(k)  
        End If
        k= k+1
     End DO

     ! restore dot
     recvBuffStart = recvBuffStart + (restrt+1)*(restrt+1)
     Do i=1,sizedot
        If(first == 1) Then
           work(recvBuffStart+i-1) = recvbuff(k)
        End If
        k = k +1
     End Do

     !restore ycurrent
    recvBuffStart = recvBuffStart+ sizedot
    Do i=1,failediter+1
       If(first == 1)Then
          work(recvBuffStart + i - 1) = recvbuff(k)
       End If
       k = k+1
    End DO

    !-------------------------------------------------------------------------
    ! [3.3] restore distributed data Xcurrent
    !-------------------------------------------------------------------------
     
     recvBuffStart = recvBuffStart + restrt 
     recvBuffend = recvBuffStart + bound_ndof-1
     Do j=ptr_Index_Intrfc(neighb),&
          ptr_Index_Intrfc(neighb+1)-1
        work(index_intrfc(j)+recvBuffStart-1) = recvbuff(k)
        k = k+1
     End Do
     !-------------------------------------------------------------------------    
     ! [3.4] restore local shared data rotsin, rotcos
     !-------------------------------------------------------------------------

     recvBuffStart = recvBuffend 
     Do i = 1,2 
        DO j =1,failediter+1
           ij = j + recvBuffStart+ (i-1)*restrt
           If (first ==1)Then
              work(ij) = recvbuff(k)
           End If
           k = k+1
        End Do
     End Do
     first = 0
  End Do
  info = 0
  
  
  !----------------------------------------------------------------------
  ! [4] Exit routine
  !----------------------------------------------------------------------
9999 Continue
  If (Associated(recvbuff)) Deallocate (recvbuff)
    If (Associated(mpiReq)) Deallocate (mpiReq)
    !Free memory
  End SUBROUTINE  XMPH_Recv_Neighbors
  

  SUBROUTINE   XMPH_FSend_Neighbors( & ! intents
       mphs,                 & ! inout
       lwork,                & ! in
       work,                 & ! in
       FailedNeighbs,        & ! in
       NbOfFailedNeighb,     & ! in
       fgmicntl,             & ! in
       failediter,           & ! in
       info                  & ! out
       )       
       
    !* Module(s) *!
    Use XMPH_maphys_type
    Use mph_error_mod        
    Use MPH_maphys_enum
    Use MPH_time_mod
    Implicit None
    Include 'mpif.h'
    
    !* Arguments *!
    Type(XMPH_maphys_t)    , Intent(inout) :: mphs
    Integer                , Intent(in )   :: lwork
    XMPH_FLOAT                             :: work(*)
    Integer                ,Intent(in  )   :: FailedNeighbs(3)
    Integer                ,Intent(in  )   :: NbOfFailedNeighb
    Integer                ,Intent(in )    :: fgmicntl(*)
    Integer                ,Intent(in )    :: failediter
    Integer                , Intent(out )  :: info
    
    !* Local variables *!

    ! Scalars
    Integer, Parameter :: PresetMPITagg = 77
    Integer :: tagg, myrank 
    Integer :: nNeighb 
    Integer :: neighb  
    Integer :: indexFailed_Neighbor
    Integer :: mpiReq
    Integer :: restrt ! restart parameter
    Integer :: bound_ndof  ! ndof of the interface
    Integer :: intrf_ndof
    Integer :: dotsize
    Integer :: failedproc
    Integer :: failedrank
    Integer :: Copy_failediter
    MPH_INT :: maxSizeIntrf
    MPH_INT ::  sendBuffSize
    MPH_INT ::  sendBuffStart
    MPH_INT ::  sendBuffend
    MPH_INT :: i,j, ij ,k 
    ! Arrays 
    Integer     , Pointer :: indexVi       (:)
    Integer     , Pointer :: ptr_Index_Intrfc (:)
    Integer     , Pointer :: Index_Intrfc(:)
    XMPH_FLOAT  , Pointer :: sendbuff             (:)
    Integer :: mpistatus(MPI_STATUS_SIZE) 
   
    !-------------------------------------------------------------------
    ![1] Initialize local arguments
    !-------------------------------------------------------------------
    ! restart parameter
    restrt  = mphs%ikeep(IKEEP_ITS_Restart) 
    ! communication  
    info = 0
    myrank= mphs%ikeep(IKEEP_MPIRANK)
    tagg  = PresetMPITagg
    !save value of faileiter
    Copy_failediter = failediter
    sendBuffSize = 0    
    ! Domain related
    nNeighb          = mphs%lc_domain%mynbvi
    indexVi          => mphs%lc_domain%myindexVi
    ptr_Index_Intrfc => mphs%lc_domain%myptrindexVi 
    indexFailed_Neighbor = -1
    bound_ndof    = mphs%lc_domain%myndofintrf
    Index_Intrfc     =>mphs%lc_domain%myindexIntrf
    intrf_ndof    = mphs%lc_domain%gballintrf    
    !! check restrt
    If (restrt .gt. intrf_ndof) restrt = intrf_ndof
    If (failediter .gt. restrt) Copy_failediter = failediter -restrt
    ! the workspace should be large enough to store the m dot-products 
    If ((fgmicntl(4).Eq.2).Or.(fgmicntl(4).Eq.3)) Then
       dotsize = restrt
    Else
       dotsize =1
    End If
    
    ! Set the ratio of lost data 
    mphs%rinfo (RINFO_RATIO_Lost) = REAL(bound_ndof*100, kind(0.d0))/intrf_ndof  

    failedrank = FailedNeighbs(1)
    !-------------------------------------------------------------------------
    ! [2] Determination of the interface with the failed neighbor 
    !-------------------------------------------------------------------------
    ! search the local index of
    ! the failed neighbor 
       Do neighb=1, nNeighb
          If (failedrank .Eq. indexVi(neighb)) Then
             indexFailed_Neighbor = neighb
             Exit  
          End If
       End Do
       
       ! k is the size of the interface, 
       k = ptr_Index_Intrfc(indexFailed_Neighbor+1) &
            - ptr_Index_Intrfc(indexFailed_Neighbor)
       
       ! compute size of the send buffer
       sendBuffSize = (restrt+1)*(Copy_failediter+1) + 2*(Copy_failediter+1)*k&
            + 5*k + 3*restrt +dotsize 
       ! Perform allocation
       Allocate( sendbuff(sendBuffSize), STAT = info)
       CHCKASSRT( info == 0, info )
       If( info < 0 ) Goto 9999
       
       !-------------------------------------------------------------------------
       ! [3] Copy and send data to the failed processor 
       !-------------------------------------------------------------------------
       
       !-------------------------------------------------------------------------
       ! [3.1] Copy distributed data (x,b, dot, r, w) 
       !------------------------------------------------------------------------
       
       k = 1
       !x,b
       DO i= 1,2
          Do j=ptr_Index_Intrfc(indexFailed_Neighbor),&
               ptr_Index_Intrfc(indexFailed_Neighbor+1)-1
             ij = index_intrfc(j) + bound_ndof*(i-1)
             sendbuff(k) = work(ij)
             k = k+1
          End Do
       End Do
       !dot
       Do i=2*bound_ndof+1,2*bound_ndof+dotsize
          work(i) = sendbuff(k)
          k = k+1
       End Do
       
       !r,w
       DO i= 3,4
          Do j=ptr_Index_Intrfc(indexFailed_Neighbor),&
               ptr_Index_Intrfc(indexFailed_Neighbor+1)-1
             ij = index_intrfc(j)+ dotsize + bound_ndof*(i-1)
             sendbuff(k) = work(ij)
             k = k+1
          End Do
       End Do
       !-------------------------------------------------------------------------   
       ! [3.2] Copy local shared data (H,ycurrent)
       !-------------------------------------------------------------------------
       
       sendBuffStart = 4*bound_ndof +dotsize +1
       
       !copy Copy_failediter+1 vectors
       Do i=sendBuffStart, sendBuffStart+& 
            (restrt+1)*(Copy_failediter+1)-1
          sendbuff(k) = work(i)
          k= k+1
       End Do
       !copy H(1,restrt+1) and ycurrent
       sendBuffStart =sendBuffStart +(restrt+1)*(restrt)
       Do i= sendBuffStart, sendBuffStart + 2*restrt
          sendbuff(k) = work(i)
          k= k+1
       End DO
       
       !-------------------------------------------------------------------------
       ! [3.3] copy distributed data Xcurrent
       !------------------------------------------------------------------------
       sendBuffStart = sendBuffStart + 2*restrt +1
       sendBuffend = sendBuffStart + bound_ndof-1
       Do j=ptr_Index_Intrfc(indexFailed_Neighbor),&
            ptr_Index_Intrfc(indexFailed_Neighbor+1)-1
          sendbuff(k) = work(index_intrfc(j)+sendBuffStart-1)
          k = k+1
       End Do
       !-------------------------------------------------------------------------    
       ! [3.4] copy local shared data rotsin, rotcos
       !-------------------------------------------------------------------------
       sendBuffStart = sendBuffend 
       Do i = 1,2 
          DO j =1,Copy_failediter+1
             ij = j + sendBuffStart+ (i-1)*restrt
          sendbuff(k) = work(ij)
          k = k+1
       End Do
    End Do
    
    !-------------------------------------------------------------------------    
    ! [3.5] copy distributed data V
    !-------------------------------------------------------------------------
    sendBuffStart=sendBuffStart+2*restrt
    DO i= 1,(Copy_failediter+1)
       Do j=ptr_Index_Intrfc(indexFailed_Neighbor),&
            ptr_Index_Intrfc(indexFailed_Neighbor+1)-1
          ij =  index_intrfc(j) + sendBuffStart + bound_ndof*(i-1)
          sendbuff(k) = work(ij)
          k = k+1
       End Do
    End Do
    
    !-------------------------------------------------------------------------    
    ! [3.6] copy distributed data Z
    !-------------------------------------------------------------------------
    
    sendBuffStart = lwork - (Copy_failediter+1)*bound_ndof
    DO i= 1,(Copy_failediter+1)
       Do j=ptr_Index_Intrfc(indexFailed_Neighbor),&
            ptr_Index_Intrfc(indexFailed_Neighbor+1)-1
          ij =  index_intrfc(j) + sendBuffStart + bound_ndof*(i-1)
          sendbuff(k) = work(ij)
          k = k+1
       End Do
    End Do
    
    !-------------------------------------------------------------------------
    ! [4 ] send data
    !-------------------------------------------------------------------------
    CALL MPI_send(                    &
         sendbuff(1),sendBuffSize,XMPH_FLOATMPI,failedrank,&
         tagg,mphs%comm,info)
    ASSRT( info == MPI_SUCCESS )
    
    info = 0
 !----------------------------------------------------------------------
 ! [5] Exit routine
 !----------------------------------------------------------------------
 
9999 Continue
 !Free memory
 If (Associated(sendbuff)) Deallocate (sendbuff)
 
End SUBROUTINE XMPH_FSend_Neighbors


SUBROUTINE XMPH_FRecv_Neighbors(& ! intents
     mphs,                 & ! in
     lwork,                & ! in
     work,                 & ! inout
     NbOfFaultyProc,       & ! in
     fgmicntl,             & ! in
     failediter,           & ! in
     info                  & ! out
     )
  
  !* Module(s) *!
  Use XMPH_maphys_type
  Use mph_error_mod        
  Use MPH_maphys_enum
  Use MPH_time_mod
  Implicit None
  Include 'mpif.h'
  
  !* Arguments *!
  Type(XMPH_maphys_t)   , Intent(inout) :: mphs
  Integer               , Intent(in )   :: lwork
  XMPH_FLOAT            , Intent(inout) :: work(*)
  Integer               , Intent(in  )  :: NbOfFaultyProc
  Integer               ,Intent(in )    :: fgmicntl(*)
  Integer               ,Intent(in )    :: failediter
  Integer               , Intent(out )  :: info
  
  !* Local variables *!
  ! Scalars
  Integer, Parameter :: PresetMPITagg = 77
  Integer :: tagg, myrank 
  Integer :: nNeighb 
  Integer :: neighb  
  Integer :: restrt ! restart parameter
  Integer :: bound_ndof  ! ndof of the interface
  Integer :: intrf_ndof
  Integer :: first
  Integer :: dotsize
  Integer :: Copy_failediter
  MPH_INT :: maxBuffSize
  MPH_INT :: maxInterSize
  MPH_INT ::  recvBuffSize
  MPH_INT ::  recvBuffStart
  MPH_INT ::  recvBuffend
  MPH_INT :: i,j, ij ,k,l 
  ! Arrays 
  Integer     , Pointer :: indexVi       (:)
  Integer     , Pointer :: ptr_Index_Intrfc (:)
  Integer     , Pointer :: Index_Intrfc     (:)
  XMPH_FLOAT  , Pointer :: recvbuff          (:)
  Integer     , Pointer :: mpiReq (:)
  Integer :: mpistatus(MPI_STATUS_SIZE) 
  


  !-------------------------------------------------------------------
  ![1] Initialize local arguments
  !-------------------------------------------------------------------
  ! restart parameter
  restrt  = mphs%ikeep(IKEEP_ITS_Restart) 
  ! communication  
  info = 0
  !save failediter 
  Copy_failediter = failediter
  myrank= mphs%ikeep(IKEEP_MPIRANK)
  tagg  = PresetMPITagg
  recvBuffSize = 0    
  ! Domain related
  nNeighb          = mphs%lc_domain%mynbvi
  indexVi          => mphs%lc_domain%myindexVi
  ptr_Index_Intrfc => mphs%lc_domain%myptrindexVi 
  bound_ndof    = mphs%lc_domain%myndofintrf
  Index_Intrfc     => mphs%lc_domain%myindexIntrf
  intrf_ndof    = mphs%lc_domain%gballintrf
  maxInterSize = 0
  first = 1
  Allocate ( mpiReq(nNeighb), STAT = info )
  CHCKASSRT( info == 0, info )
  If( info < 0 ) Goto 9999
  !! check restrt
  If (restrt .ge. intrf_ndof)restrt = intrf_ndof
  If (failediter .gt. restrt)Copy_failediter = failediter - restrt
  
  ! Set the ratio of lost data 
  mphs%rinfo (RINFO_RATIO_Lost) = REAL(bound_ndof*100, kind(0.d0))/intrf_ndof  
  ! set number of processors involved in recovery 
  mphs%rinfo (RINFO_PROC_INVOLVED) = REAL(nNeighb+1, kind(0.d0))

  ! the workspace should be large enough to store the m dot-products 
  If ((fgmicntl(4).Eq.2).Or.(fgmicntl(4).Eq.3)) Then
     dotsize = restrt
  Else
     dotsize =1
  End If
  !-------------------------------------------------------------------------
  ! [2] Detect message from neighbor and allocate receive buffer 
  !-------------------------------------------------------------------------

    !search maximum of interfaces size
     Do neighb=1,nNeighb 
        k = ptr_Index_Intrfc(neighb+1)-ptr_Index_Intrfc(neighb)
        If (k > maxInterSize) Then
           maxInterSize = k
        End If
     End Do
        maxBuffSize = (restrt+1)*(Copy_failediter+1) + 2*(Copy_failediter+1)*maxInterSize&
             + 5*maxInterSize +3*restrt +dotsize

        ! Perform allocation
        Allocate (recvbuff(maxBuffSize),STAT = info )
        CHCKASSRT( info == 0, info )
        If( info < 0 ) Goto 9999

  Do l=1,nNeighb 
     ! Get a message
     Call MPI_Probe(MPI_ANY_SOURCE,tagg,mphs%comm,mpiStatus,info)
     ASSRT( info == MPI_SUCCESS )
     
     Call MPI_Get_count(mpiStatus,XMPH_FLOATMPI,recvBuffsize,info)
     ASSRT( info == MPI_SUCCESS )

     ! Detect which neighbor this message corresponds to
     Do neighb=1,nNeighb
	If ( indexVi(neighb)  ==  mpiStatus(MPI_SOURCE) ) Exit 
     End Do

     !-------------------------------------------------------------------------
     ! [3] receive and restore data 
     !-------------------------------------------------------------------------

     Call MPI_recv(recvbuff(1),recvBuffSize,&
          XMPH_FLOATMPI,mpiStatus(MPI_SOURCE),&
          tagg,mphs%comm,MPI_STATUS_IGNORE,info)
     ASSRT( info == MPI_SUCCESS )
     !-------------------------------------------------------------------------
     ! [3.1] restore distributed data (x,b, dot, r, w) 
     !-------------------------------------------------------------------------
     !x, b
     k = 1         
     Do i=1,2
        Do j=ptr_Index_Intrfc(neighb),&
             ptr_Index_Intrfc(neighb+1)-1
           ij = index_intrfc(j) + bound_ndof*(i-1)
           work(ij) = recvbuff(k)
           k = k +1
        End Do
     End Do
     !dot
     Do i = 2*bound_ndof+1, 2*bound_ndof+dotsize
        work(i) = recvbuff(k)
        k = k+1    
     End Do
     !r, w
     Do i=3,4
        Do j=ptr_Index_Intrfc(neighb),&
             ptr_Index_Intrfc(neighb+1)-1
           ij = dotsize+ index_intrfc(j) + bound_ndof*(i-1)
           work(ij) = recvbuff(k)
           k = k +1
        End Do
     End Do     

     !-------------------------------------------------------------------------   
     ! [3.2] restore local shared data (H,ycurrent)  
     !-------------------------------------------------------------------------
     
    recvBuffStart = 4*bound_ndof +1+dotsize

     ! restore (Copy_failediter+1) vectors 
     Do i=recvBuffStart, recvBuffStart+&
          (restrt+1)*(Copy_failediter+1)-1
        if(first == 1) Then   
           work(i) =   recvbuff(k)  
        End If
        k = k+1
     End Do
     ! restore H(1,m+1) and ycurrent
    recvBuffStart =recvBuffStart +(restrt+1)*(restrt)
    Do i= recvBuffStart, recvBuffStart + 2*restrt
        if(first == 1) Then   
           work(i) =   recvbuff(k)  
        End If
        k= k+1
     End DO

    !-------------------------------------------------------------------------
    ! [3.3] restore distributed data Xcurrent
    !-------------------------------------------------------------------------
     
     recvBuffStart = recvBuffStart + 2*restrt +1 
     recvBuffend = recvBuffStart + bound_ndof-1
     Do j=ptr_Index_Intrfc(neighb),&
          ptr_Index_Intrfc(neighb+1)-1
        work(index_intrfc(j)+recvBuffStart-1) = recvbuff(k)
        k = k+1
     End Do
     !-------------------------------------------------------------------------    
     ! [3.4] restore local shared data rotsin, rotcos
     !-------------------------------------------------------------------------
     recvBuffStart = recvBuffend 
     Do i = 1,2 
        DO j =1,Copy_failediter+1
           ij = j + recvBuffStart+ (i-1)*restrt
           If (first ==1)Then
              work(ij) = recvbuff(k)
           End If
           k = k+1
        End Do
     End Do

     !-------------------------------------------------------------------------    
     ! [3.4] restore distributed data V
     !-------------------------------------------------------------------------
     recvbuffstart = recvbuffstart +2*restrt
     DO i= 1,(Copy_failediter+1)
        Do j=ptr_Index_Intrfc(neighb),&
             ptr_Index_Intrfc(neighb+1)-1
           ij =  index_intrfc(j) + recvBuffStart + bound_ndof*(i-1)
           work(ij) = recvbuff(k)
           k = k+1
        End Do
     End Do
     
     !-------------------------------------------------------------------------    
     ! [3.6] restore distributed data Z
     !-------------------------------------------------------------------------
    recvBuffStart = lwork - (Copy_failediter+1)*bound_ndof
    DO i= 1,(Copy_failediter+1)
       Do j=ptr_Index_Intrfc(neighb),&
            ptr_Index_Intrfc(neighb+1)-1
          ij =  index_intrfc(j) + recvBuffStart + bound_ndof*(i-1)
          work(ij) = recvbuff(k)
          k = k+1
       End Do
    End Do
     first = 0
  End Do
  info = 0
  !----------------------------------------------------------------------
  ! [4] Exit routine
  !----------------------------------------------------------------------
9999 Continue
    If (Associated(recvbuff)) Deallocate (recvbuff)
    If (Associated(mpiReq)) Deallocate (mpiReq)
  !Free memory
End SUBROUTINE  XMPH_FRecv_Neighbors



SUBROUTINE  XMPH_ExchContrib(& ! Intent 
     mphs,                      & ! in
     work,                      & ! inout
     FailedNeighbs,             & ! in
     NbOfFailedNeighb,          & ! in
     Iamfailed                  & ! in 
     )
  
  !* Module(s) *!
  Use XMPH_maphys_type
  Use mph_error_mod        
  Use MPH_maphys_enum
  Use XMPH_dense_matrix_mod  
  Implicit None
  Include 'mpif.h'
  
  !* Arguments *!
  Type(XMPH_maphys_t)   , Intent(inout ) :: mphs
  XMPH_FLOAT            , Intent(inout ) :: work(*) 
  Integer               , Intent(in )  :: FailedNeighbs(3)
  Integer               , Intent(in )  :: NbOfFailedNeighb
  Integer               , Intent(in )  :: Iamfailed

    ! Scalars
    Integer, Parameter :: PresetMPITagg = 77
    Integer :: tagg, myrank, tagg2 
    Integer :: nNeighb
    Integer :: neighb  
    Integer :: indexFailed_Neighbor
    Integer :: bound_ndof  ! ndof of the interface
    Integer :: intrf_ndof
    Integer :: failedrank
    Integer :: failedproc
    Integer :: info
    Integer :: bufferSize
    MPH_INT :: maxSizeIntrf
    MPH_INT ::  sendBuffSize
    MPH_INT ::  sendBuffStart
    MPH_INT ::  sendBuffend
    MPH_INT ::  recvBuffSize
    MPH_INT ::  maxBuffSize
    MPH_INT :: i,j, ij ,k ,l
    ! Arrays 
    Integer     , Pointer :: indexVi       (:)
    Integer     , Pointer :: ptr_Index_Intrfc (:)
    Integer     , Pointer :: Index_Intrfc(:)
    Integer     , Pointer :: mpiReq(:)
    XMPH_FLOAT  , Pointer :: sendbuff          (:)
    XMPH_FLOAT  , Pointer :: recvbuff      (:)
    Integer :: mpistatus(MPI_STATUS_SIZE) 
    ! Derived types
    Type(XMPH_dense_matrix_t) :: x
    Type(XMPH_dense_matrix_t) :: z

    !-------------------------------------------------------------------
    ![1] Initialize local arguments
    !-------------------------------------------------------------------
    ! communication  
    info = 0
    myrank= mphs%ikeep(IKEEP_MPIRANK)
    tagg  = PresetMPITagg
    tagg2 = 33
    sendBuffSize = 0    
    ! Domain related
    nNeighb          = mphs%lc_domain%mynbvi
    indexVi          => mphs%lc_domain%myindexVi
    ptr_Index_Intrfc => mphs%lc_domain%myptrindexVi 
    indexFailed_Neighbor = -1
    bound_ndof    = mphs%lc_domain%myndofintrf
    Index_Intrfc     =>mphs%lc_domain%myindexIntrf
    intrf_ndof    = mphs%lc_domain%gballintrf    
    maxSizeIntrf = 0    
    Sendbuffstart = 1

    Allocate ( mpiReq(10), STAT = info )
    CHCKASSRT( info == 0, info )
    If( info < 0 ) Goto 9999
    
    Call XMPH_DM_Create(x, bound_ndof, 1, bound_ndof, info)
    MPH_ONFAILURE_GOTO9999(info)
    Call XMPH_DM_Create(z, bound_ndof, 1, bound_ndof, info)
    MPH_ONFAILURE_GOTO9999(info)
    
    !search maximum of interfaces size
    maxBuffSize = 2*bound_ndof 

    ! set ratio of lost data
    mphs%rinfo (RINFO_RATIO_Lost)= REAL(bound_ndof*100, kind(0.d0))/intrf_ndof  
    ! Perform allocation
    Allocate (recvbuff(maxBuffSize), STAT = info )
    CHCKASSRT( info == 0, info )
    If( info < 0 ) Goto 9999
    Allocate (sendbuff(maxBuffSize), STAT = info )
    CHCKASSRT( info == 0, info )
    If( info < 0 ) Goto 9999

    !-------------------------------------------------------------------------
    ! [2 ] Neighbors of failed processors prepare and send contribution
    !-------------------------------------------------------------------------
    If(Iamfailed .Eq. 0 ) Then              
       Do failedproc=1,NbOfFailedNeighb
          failedrank=FailedNeighbs(failedproc)
          ! search the local index of the failed neighbor 
          Do neighb=1, nNeighb
             If (failedrank .Eq. indexVi(neighb)) Then
                indexFailed_Neighbor = neighb
                Exit  
             End If
          End Do
          
          ! k is the size of the interface, 
          k = ptr_Index_Intrfc(indexFailed_Neighbor+1) &
               - ptr_Index_Intrfc(indexFailed_Neighbor)       
          sendBuffSize = 2*k
          !Copy the current solution to x_local
          Call XMPH_ARITHcopy &
               (bound_ndof, work(1),1, x%v(1), 1)

          !Set entries shared with failed proc to zero       
          Do i=ptr_Index_Intrfc(indexFailed_Neighbor),&
               ptr_Index_Intrfc(indexFailed_Neighbor+1)-1
             x%v(index_intrfc(i)) = XMPH_FLOATZERO
          End Do
          
          ! Schur * x Product 
          Call XMPH_DM_VectorProduct(mphs%dm_schur, x, z, info )
          MPH_ONFAILURE_GOTO9999(info)
          
          !copy x and rhs contribution to sendbuffer    
          K= Sendbuffstart
          Do i=ptr_Index_Intrfc(indexFailed_Neighbor),&
               ptr_Index_Intrfc(indexFailed_Neighbor+1)-1
             sendbuff(k) = work(index_intrfc(i))
             k= k+1
          End Do
          
          Do i=ptr_Index_Intrfc(indexFailed_Neighbor),&
            ptr_Index_Intrfc(indexFailed_Neighbor+1)-1
             sendbuff(k) = z%v(index_intrfc(i))
             k=k+1
          End Do
          
          !-------------------------------------------------------------------------
          ! [2.1 ] send data
          !-------------------------------------------------------------------------
          CALL MPI_Isend(                    &
               sendbuff(Sendbuffstart),sendBuffSize,XMPH_FLOATMPI,failedrank,&
               tagg,mphs%comm, mpiReq(failedproc), info)
          ASSRT( info == MPI_SUCCESS )
          Sendbuffstart = k
       End Do
       Call MPI_Waitall(NbOfFailedNeighb, mpiReq, MPI_STATUSES_IGNORE,info)
       ASSRT( info == MPI_SUCCESS )
    End IF

    !-------------------------------------------------------------------------
    ! [3] Receive contribution from non failed neighbors
    !-------------------------------------------------------------------------

    If(( Iamfailed .Eq. 1) .And. ((nNeighb-NbOfFailedNeighb) .Gt. 0) ) Then
       !Copy rhs to 2*bound_ndof+1
       Do i=1,bound_ndof
          work(2*bound_ndof+i) = work(bound_ndof+i)
       End DO

       Do l=1,nNeighb-NbOfFailedNeighb
          ! Get a message
          Call MPI_Probe(MPI_ANY_SOURCE,tagg,mphs%comm,mpiStatus,info)
          ASSRT( info == MPI_SUCCESS )
          
          Call MPI_Get_count(mpiStatus,XMPH_FLOATMPI,recvBuffsize,info)
          ASSRT( info == MPI_SUCCESS )
          
          ! Detect which neighbor this message corresponds to
          Do neighb=1,nNeighb
             If ( indexVi(neighb)  ==  mpiStatus(MPI_SOURCE) ) Exit 
          End Do
          
          !-------------------------------------------------------------------------
          ! [3.1] receive and restore data 
          !-------------------------------------------------------------------------

          Call MPI_recv(recvbuff(1),recvBuffSize,&
               XMPH_FLOATMPI,mpiStatus(MPI_SOURCE),&
               tagg,mphs%comm,MPI_STATUS_IGNORE,info)
          ASSRT( info == MPI_SUCCESS )

          !restore x
          k=1
          Do i=ptr_Index_Intrfc(neighb),&
               ptr_Index_Intrfc(neighb+1)-1     
             work(index_intrfc(i)) = recvbuff(k)
             k=k+1
          End Do
          ! update rhs
          Do i=ptr_Index_Intrfc(neighb),&
               ptr_Index_Intrfc(neighb+1)-1     
             j = 2*bound_ndof+index_intrfc(i)
             work(j) = work(j) - recvbuff(k)
             k=k+1
          End Do
        End Do
     End If

    !--------------------------------------------------------           
    ! [4] data exchange between failed procs
    !---------------------------------------------------------
    
    If((Iamfailed .Eq. 1) .AND. (NbOfFailedNeighb .Gt. 0) )Then
       failedrank=FailedNeighbs(NbOfFailedNeighb)
          
       ! search the local index of the failed neighbor 
          Do neighb=1, nNeighb
             If (failedrank .Eq. indexVi(neighb)) Then
                indexFailed_Neighbor = neighb
                Exit  
             End If
          End Do
          ! k is the size of the interface, 
          bufferSize = ptr_Index_Intrfc(indexFailed_Neighbor+1) &
               - ptr_Index_Intrfc(indexFailed_Neighbor)       
          
          !Copy the current solution to x_local
          Call XMPH_ARITHcopy &
            (bound_ndof, work(1),1, x%v(1), 1)
          !Set entries shared with failed proc to zero       
          Do i=ptr_Index_Intrfc(indexFailed_Neighbor),&
               ptr_Index_Intrfc(indexFailed_Neighbor+1)-1
             x%v(index_intrfc(i)) = XMPH_FLOATZERO
          End Do
          
          ! Schur * x Product 
          Call XMPH_DM_VectorProduct(mphs%dm_schur, x, z, info )
          MPH_ONFAILURE_GOTO9999(info)

          !copy rhs contribution    
          k=1
          Do i=ptr_Index_Intrfc(indexFailed_Neighbor),&
               ptr_Index_Intrfc(indexFailed_Neighbor+1)-1
             sendbuff(k) = z%v(index_intrfc(i))
             k=k+1
          End Do

          !send
         
          CALL MPI_Isend(                    &
               sendbuff(1),bufferSize,XMPH_FLOATMPI,failedrank,&
               tagg2,mphs%comm,mpiReq(1), info)
          ASSRT( info == MPI_SUCCESS )
          
          Call MPI_Irecv(recvbuff(1),bufferSize,&
               XMPH_FLOATMPI,failedrank,&
               tagg2,mphs%comm,mpiReq(2) ,info)
          ASSRT( info == MPI_SUCCESS )              
          Call MPI_Waitall(2, mpiReq, MPI_STATUSES_IGNORE, info)
          ASSRT( info == MPI_SUCCESS )              
          k=1  
          Do i=ptr_Index_Intrfc(indexFailed_Neighbor),&
               ptr_Index_Intrfc(indexFailed_Neighbor+1)-1     
             j = 2*bound_ndof+index_intrfc(i)
             work(j) = work(j) - recvbuff(k)
             k=k+1
          End Do
       End If
9999 Continue    

  If (Associated(x%v)) Deallocate  (x%v)
  If (Associated(z%v)) Deallocate  (z%v)
  If (Associated(recvbuff)) Deallocate (recvbuff)
  If (Associated(mpiReq)) Deallocate  (mpiReq)
  If (Associated(sendbuff)) Deallocate (sendbuff)

End SUBROUTINE XMPH_ExchContrib


SUBROUTINE  XMPH_RestoreX(& ! Intent 
     mphs,                      & ! in
     work,                      & ! inout
     FailedNeighbs,             & ! in
     NbOfFailedNeighb,          & ! in
     Iamfailed                  & ! in 
     )
  
  !* Module(s) *!
  Use XMPH_maphys_type
  Use mph_error_mod        
  Use MPH_maphys_enum
  Use XMPH_dense_matrix_mod  
  Implicit None
  Include 'mpif.h'
  
  !* Arguments *!
  Type(XMPH_maphys_t)   , Intent(inout ) :: mphs
  XMPH_FLOAT            , Intent(inout ) :: work(*) 
  Integer               , Intent(in )  :: FailedNeighbs(3)
  Integer               , Intent(in )  :: NbOfFailedNeighb
  Integer               , Intent(in )  :: Iamfailed

    ! Scalars
    Integer, Parameter :: PresetMPITagg = 77
    Integer :: tagg, myrank, tagg2 
    Integer :: nNeighb
    Integer :: neighb  
    Integer :: indexFailed_Neighbor
    Integer :: bound_ndof  ! ndof of the interface
    Integer :: intrf_ndof
    Integer :: failedrank
    Integer :: failedproc
    Integer :: info
    Integer :: bufferSize
    MPH_INT :: maxSizeIntrf
    MPH_INT ::  sendBuffSize
    MPH_INT ::  sendBuffStart
    MPH_INT ::  recvBuffSize
    MPH_INT ::  maxBuffSize
    MPH_INT :: i,j, ij ,k ,l
    ! Arrays 
    Integer     , Pointer :: indexVi       (:)
    Integer     , Pointer :: ptr_Index_Intrfc (:)
    Integer     , Pointer :: Index_Intrfc(:)
    Integer     , Pointer :: mpiReq(:)
    XMPH_FLOAT  , Pointer :: sendbuff          (:)
    XMPH_FLOAT  , Pointer :: recvbuff      (:)
    Integer :: mpistatus(MPI_STATUS_SIZE) 

    !-------------------------------------------------------------------
    ![1] Initialize local arguments
    !-------------------------------------------------------------------
    ! communication  
    info = 0
    myrank= mphs%ikeep(IKEEP_MPIRANK)
    tagg  = PresetMPITagg
    tagg2 = 33
    sendBuffSize = 0    
    ! Domain related
    nNeighb          = mphs%lc_domain%mynbvi
    indexVi          => mphs%lc_domain%myindexVi
    ptr_Index_Intrfc => mphs%lc_domain%myptrindexVi 
    indexFailed_Neighbor = -1
    bound_ndof    = mphs%lc_domain%myndofintrf
    Index_Intrfc     =>mphs%lc_domain%myindexIntrf
    intrf_ndof    = mphs%lc_domain%gballintrf    
    maxSizeIntrf = 0    
    Sendbuffstart = 1

    Allocate ( mpiReq(10), STAT = info )
    CHCKASSRT( info == 0, info )
    If( info < 0 ) Goto 9999
    
    !search maximum of interfaces size
    maxBuffSize = bound_ndof 

    ! Perform allocation
    Allocate (recvbuff(maxBuffSize), STAT = info )
    CHCKASSRT( info == 0, info )
    If( info < 0 ) Goto 9999
    Allocate (sendbuff(maxBuffSize), STAT = info )
    CHCKASSRT( info == 0, info )
    If( info < 0 ) Goto 9999

    ! set ratio of lost data
    mphs%rinfo (RINFO_RATIO_Lost)= REAL(bound_ndof*100, kind(0.d0))/intrf_ndof  

    !-------------------------------------------------------------------------
    ! [2 ] Neighbors of failed processors prepare and send contribution
    !-------------------------------------------------------------------------
    If(Iamfailed .Eq. 0 ) Then              
       Do failedproc=1,NbOfFailedNeighb
          failedrank=FailedNeighbs(failedproc)
          ! search the local index of the failed neighbor 
          Do neighb=1, nNeighb
             If (failedrank .Eq. indexVi(neighb)) Then
                indexFailed_Neighbor = neighb
                Exit  
             End If
          End Do
          
          ! k is the size of the interface, 
          k = ptr_Index_Intrfc(indexFailed_Neighbor+1) &
               - ptr_Index_Intrfc(indexFailed_Neighbor)       
          sendBuffSize = k

          !copy x 
          K= Sendbuffstart
          Do i=ptr_Index_Intrfc(indexFailed_Neighbor),&
               ptr_Index_Intrfc(indexFailed_Neighbor+1)-1
             sendbuff(k) = work(index_intrfc(i))
             k= k+1
          End Do
          
          !-------------------------------------------------------------------------
          ! [2.1 ] send data
          !-------------------------------------------------------------------------
          CALL MPI_Isend(                    &
               sendbuff(Sendbuffstart),sendBuffSize,XMPH_FLOATMPI,failedrank,&
               tagg,mphs%comm, mpiReq(failedproc), info)
          ASSRT( info == MPI_SUCCESS )
          Sendbuffstart = k
       End Do
       Call MPI_Waitall(NbOfFailedNeighb, mpiReq, MPI_STATUSES_IGNORE,info)
       ASSRT( info == MPI_SUCCESS )
    End IF

    !-------------------------------------------------------------------------
    ! [3] Receive contribution from non failed neighbors
    !-------------------------------------------------------------------------

    If(( Iamfailed .Eq. 1) .And. ((nNeighb-NbOfFailedNeighb) .Gt. 0) ) Then

       Do l=1,nNeighb-NbOfFailedNeighb
          ! Get a message
          Call MPI_Probe(MPI_ANY_SOURCE,tagg,mphs%comm,mpiStatus,info)
          ASSRT( info == MPI_SUCCESS )
          
          Call MPI_Get_count(mpiStatus,XMPH_FLOATMPI,recvBuffsize,info)
          ASSRT( info == MPI_SUCCESS )
          
          ! Detect which neighbor this message corresponds to
          Do neighb=1,nNeighb
             If ( indexVi(neighb)  ==  mpiStatus(MPI_SOURCE) ) Exit 
          End Do
          
          !-------------------------------------------------------------------------
          ! [3.1] receive and restore data 
          !-------------------------------------------------------------------------

          Call MPI_recv(recvbuff(1),recvBuffSize,&
               XMPH_FLOATMPI,mpiStatus(MPI_SOURCE),&
               tagg,mphs%comm,MPI_STATUS_IGNORE,info)
          ASSRT( info == MPI_SUCCESS )

          !restore x
          k=1
          Do i=ptr_Index_Intrfc(neighb),&
               ptr_Index_Intrfc(neighb+1)-1     
             work(index_intrfc(i)) = recvbuff(k)
             k=k+1
          End Do
        End Do
     End If
9999 Continue    

  If (Associated(recvbuff)) Deallocate (recvbuff)
  If (Associated(mpiReq)) Deallocate  (mpiReq)
  If (Associated(sendbuff)) Deallocate (sendbuff)

End SUBROUTINE XMPH_RestoreX






SUBROUTINE XMPH_RecvRhsContrib(& !Intents
     mphs,                     & ! in
     work,                     & ! in
     NbOfFaultyProc            & ! in
     )

  Use XMPH_maphys_type
  Use mph_error_mod        
  Use MPH_maphys_enum
  Implicit None
  Include 'mpif.h'
  
  !* Arguments *!
  Type(XMPH_maphys_t)   , Intent(in ) :: mphs
  XMPH_FLOAT            , Intent(inout ) :: work(*) 
  Integer               , Intent(in )  :: NbOfFaultyProc

    ! Scalars
    Integer, Parameter :: PresetMPITagg = 77
    Integer :: tagg, myrank 
    Integer :: nNeighb 
    Integer :: neighb  
    Integer :: indexFailed_Neighbor
    Integer :: bound_ndof  ! ndof of the interface
    Integer :: intrf_ndof
    Integer :: failedrank
    Integer :: failedproc
    Integer :: info
    MPH_INT :: maxSizeIntrf
    MPH_INT ::  recvBuffSize
    MPH_INT ::  maxBuffSize
    MPH_INT ::  recvBuffStart
    MPH_INT ::  recvBuffend
    MPH_INT :: i,j, ij ,k,l 
    ! Arrays 
    Integer     , Pointer :: indexVi       (:)
    Integer     , Pointer :: mpiReq(:)
    Integer     , Pointer :: ptr_Index_Intrfc (:)
    Integer     , Pointer :: Index_Intrfc(:)
    XMPH_FLOAT  , Pointer :: recvbuff      (:)
    Integer :: mpistatus(MPI_STATUS_SIZE) 
   
    !-------------------------------------------------------------------
    ![1] Initialize local arguments
    !-------------------------------------------------------------------
    ! communication  
    info = 0
    myrank= mphs%ikeep(IKEEP_MPIRANK)
    tagg  = PresetMPITagg
    recvBuffSize = 0    
    ! Domain related
    nNeighb          = mphs%lc_domain%mynbvi
    indexVi          => mphs%lc_domain%myindexVi
    ptr_Index_Intrfc => mphs%lc_domain%myptrindexVi 
    indexFailed_Neighbor = -1
    bound_ndof    = mphs%lc_domain%myndofintrf
    Index_Intrfc     =>mphs%lc_domain%myindexIntrf
    intrf_ndof    = mphs%lc_domain%gballintrf    
    maxSizeIntrf = 0

Allocate ( mpiReq(nNeighb), STAT = info )
  CHCKASSRT( info == 0, info )
  If( info < 0 ) Goto 9999

  !-------------------------------------------------------------------------
  ! [2] Detect message from neighbor and allocate receive buffer 
  !-------------------------------------------------------------------------
    
    !search maximum of interfaces size
     Do neighb=1,nNeighb 
        k = ptr_Index_Intrfc(neighb+1)-ptr_Index_Intrfc(neighb)
        If (k > maxSizeIntrf) Then
           maxSizeIntrf = k
        End If
     End Do
     maxBuffSize = 2*maxSizeIntrf
     ! Perform allocation
     Allocate (recvbuff(maxBuffSize),STAT = info )
     CHCKASSRT( info == 0, info )
     If( info < 0 ) Goto 9999


  Do l=1,nNeighb
     ! Get a message
     Call MPI_Probe(MPI_ANY_SOURCE,tagg,mphs%comm,mpiStatus,info)
     ASSRT( info == MPI_SUCCESS )
     
     Call MPI_Get_count(mpiStatus,XMPH_FLOATMPI,recvBuffsize,info)
     ASSRT( info == MPI_SUCCESS )

     ! Detect which neighbor this message corresponds
     Do neighb=1,nNeighb
	If ( indexVi(neighb)  ==  mpiStatus(MPI_SOURCE) ) Exit 
     End Do

     !-------------------------------------------------------------------------
     ! [3] receive and restore data 
     !-------------------------------------------------------------------------

     Call MPI_recv(recvbuff(1),recvBuffSize,&
          XMPH_FLOATMPI,mpiStatus(MPI_SOURCE),&
          tagg,mphs%comm,mpiReq(l),info)
     ASSRT( info == MPI_SUCCESS )
     !--------------------------------------------------------           
     ! [3.1] restore distributed data x and rhs contribution
     !---------------------------------------------------------

     k=1
     Do i=ptr_Index_Intrfc(neighb),&
          ptr_Index_Intrfc(neighb+1)-1     
        work(index_intrfc(i)) = recvbuff(k)
        k=k+1
     End Do

     Do i=ptr_Index_Intrfc(neighb),&
          ptr_Index_Intrfc(neighb+1)-1     
        j = bound_ndof+1+index_intrfc(i)
        work(j) = work(j) - recvbuff(k)
        k=k+1
     End Do      
  End Do


  
9999 continue
  If (Associated(recvbuff)) Deallocate (recvbuff)
  If (Associated(mpiReq)) Deallocate  (mpiReq)
End SUBROUTINE XMPH_RecvRhsContrib


SUBROUTINE XMPH_INTERPOLATE( & ! intents
       mphs,                 & ! in 
       work,                 & ! in 
       ArrayOfFaultyPorc,    & ! in
       NbOfFaultyProc,       & ! in
       FailedNeighbs,        & ! in
       NbOfFailedNeighb,     & ! in
       info                  & ! out 
       )


  Use XMPH_maphys_type
  Use XMPH_dense_matrix_mod  
  Use XMPH_dds_mod
  Use XMPH_sds_mod
  Use mph_error_mod        
  Use MPH_maphys_enum
  Implicit None
  Include 'mpif.h'
  
  !* Arguments *!
  Type(XMPH_maphys_t)   , Intent(inout ) :: mphs
  XMPH_FLOAT            , Intent(inout ) :: work(*) 
  Integer               , Intent(in  )  :: ArrayOfFaultyPorc(3)
  Integer               , Intent(in  )  :: NbOfFaultyProc
  Integer               , Intent(out  )  :: FailedNeighbs(3)
  Integer               , Intent(out )  :: NbOfFailedNeighb
  Integer               , Intent(out ) :: info    

    Integer :: nNeighb 
    Integer :: neighb  
    Integer :: nb_Common_neighbor
    Integer :: Common
    Integer :: Common_neighbor(100)
    Integer :: Array_of_Common_index(100000)
    Integer :: nb_of_Common_index
    Integer :: indexCommon_neighbor
    Integer :: common_neighb
    Integer :: indexFailed_Neighbor
    Integer :: bound_ndof  ! ndof of the interface
    Integer :: intrf_ndof
    Integer :: failedrank
    Integer :: failedproc
    Integer :: mpiReq(2)
    Integer :: tagg, myrank
    Integer :: count_commoun
    Integer :: Pcd_Strategy
    MPH_INT :: i,j,k
    ! Arrays 
    Integer     , Pointer :: indexVi       (:)
    Integer     , Pointer :: ptr_Index_Intrfc (:)
    Integer     , Pointer :: Index_Intrfc(:)

    XMPH_FLOAT  , Pointer :: recvbuff      (:)
    XMPH_FLOAT  , Pointer :: sendbuff      (:)
    ! Derived types
    Type(XMPH_dense_matrix_t) :: z
    count_commoun = 0
    nb_of_Common_index = 0
    nb_Common_neighbor = 0
    info =0
    tagg = 777
    myrank= mphs%ikeep(IKEEP_MPIRANK)
    nNeighb          = mphs%lc_domain%mynbvi
    indexVi          => mphs%lc_domain%myindexVi
    ptr_Index_Intrfc => mphs%lc_domain%myptrindexVi 
    indexFailed_Neighbor = -1
    bound_ndof    = mphs%lc_domain%myndofintrf
    Index_Intrfc     =>mphs%lc_domain%myindexIntrf
    intrf_ndof    = mphs%lc_domain%gballintrf    
    Pcd_Strategy  = mphs%ikeep(IKEEP_PCD_Strategy)
    Call XMPH_DM_Create(z, bound_ndof, 1, bound_ndof, info)
    MPH_ONFAILURE_GOTO9999(info)


    
    ! Copy the rhs to z
    Call XMPH_ARITHcopy &
         (bound_ndof, work(2*bound_ndof+1), 1,z%v(1), 1)       

    
    If ( Pcd_Strategy .Eq. PCD_STRATEGY_isLocalExact ) Then
       Call XMPH_DDS_Solve(  &            
            mphs%dls_precond_schur%dds , & 
            mphs%dls_precond_schur%dm_A, &
            z, info)
       CHCKRINFO(info)
       If (info < 0) Goto 9999       
    Else
       Call XMPH_SDS_Solve_RHS(  &
            mphs%sls_precond_schur%sds, &
            z%v, z%n, z%ld,    &
            info)
       CHCKRINFO(info)
       If (info < 0) Goto 9999
    End If

    Call XMPH_IsNeighbor(mphs, ArrayOfFaultyPorc, NbOfFaultyProc,&
         FailedNeighbs, NbOfFailedNeighb,info)
    failedrank=FailedNeighbs(NbOfFailedNeighb)       

    
         !search the local index of the failed neighbor 
         Do neighb=1, nNeighb
            If (failedrank .Eq. indexVi(neighb)) Then
               indexFailed_Neighbor = neighb
               Exit  
            End If
         End Do
         
         k = ptr_Index_Intrfc(indexFailed_neighbor+1)-&
              ptr_Index_Intrfc(indexFailed_neighbor)       
         
         ! set ratio of interpolated data 
         mphs%rinfo (RINFO_RATIO_Interp)= REAL(k*100, kind(0.d0))/intrf_ndof  

         Allocate (recvbuff(k), sendbuff(k), STAT = info )
         CHCKASSRT( info == 0, info )
         If( info < 0 ) Goto 9999        
         
         !set send data into the buffer
         j=1
         Do i=ptr_Index_Intrfc(indexFailed_neighbor),&
              ptr_Index_Intrfc(indexFailed_neighbor+1)-1
            work(2*bound_ndof+ index_intrfc(i)) = 0.5*z%v(index_intrfc(i))
            sendbuff(j)  = 0.5*z%v(index_intrfc(i))
            j=j+1
         End Do

         Call MPI_Isend(                    &
              sendbuff(1),k,XMPH_FLOATMPI,failedrank,&
              tagg,mphs%comm, mpiReq(1),info)
         ASSRT( info == MPI_SUCCESS )       

         Call MPI_Irecv(recvbuff(1),k,&
              XMPH_FLOATMPI,failedrank,&
              tagg,mphs%comm,mpiReq(2),info)
         ASSRT( info == MPI_SUCCESS )
         Call MPI_Waitall(2, mpiReq,MPI_STATUSES_IGNORE, info)
         ASSRT( info == MPI_SUCCESS )

         info=0
        Call MPH_Common_Neighbor(mphs, Failedrank, Common_neighbor, nb_Common_neighbor)
         If(nb_Common_neighbor .Eq. 0) Then
            j=1
            Do i=ptr_Index_Intrfc(indexfailed_neighbor),&
                 ptr_Index_Intrfc(indexfailed_neighbor+1)-1
               work(index_intrfc(i)) =  work(2*bound_ndof+index_intrfc(i)) + recvbuff(j)
               j=j+1
            End Do
         Else
            !Get index of common_neighbor            
            Do Common=1,nb_common_neighbor
               Do neighb=1, nNeighb
                  If (Common_neighbor(Common) .Eq. indexVi(neighb)) Then
                     indexCommon_Neighbor = neighb
                     Exit  
                  End If
               End Do
               !Fill array of common index

               Do i=ptr_Index_Intrfc(indexfailed_neighbor),&
                    ptr_Index_Intrfc(indexfailed_neighbor+1)-1
                  Do k=ptr_Index_Intrfc(indexCommon_neighbor),&
                       ptr_Index_Intrfc(indexCommon_neighbor+1)-1
                     If (index_intrfc(i) .Eq. index_intrfc(k)) Then
                        Call  MPH_Not_contain(Array_of_Common_index, &
                             nb_of_Common_index, index_intrfc(i), info)
                        If(info .Eq. 1) Then
                           nb_of_Common_index = nb_of_Common_index +1
                        Array_of_Common_index(nb_of_Common_index) = index_intrfc(i)
                     End If
                     exit
                     End If
                  End Do
               End DO
            End DO
            !restore x
            j=1
            Do i=ptr_Index_Intrfc(indexfailed_neighbor),&
                 ptr_Index_Intrfc(indexfailed_neighbor+1)-1
               Call  MPH_Not_contain(Array_of_Common_index, &
                    nb_of_Common_index, index_intrfc(i), info)              
               If (info .Eq. 1) Then
                  work(index_intrfc(i)) = work(2*bound_ndof+index_intrfc(i)) + recvbuff(j)
               End If
               j=j+1
            End Do
         End If

         !!Dbg     
!           If(myrank .Eq. ArrayOfFaultyPorc(1) ) Then
!              open(unit=771,file="Fault1.txt",action="write")
!              write(771,*) "I am proc", myrank, "my neigbhors" 
!             Do neighb=1, nNeighb
!                write(771,*) indexVi(neighb)
!             End Do
!             if(nb_Common_neighbor > 0 ) Then
!                Do Common=1,nb_Common_neighbor
!                   Do neighb=1, nNeighb
!                      If (Common_neighbor(Common) .Eq. indexVi(neighb)) Then
!                         indexCommon_Neighbor = neighb
!                         Exit  
!                   End If
!                End Do
               
!                write(771,*) "common is", Common_neighbor(Common)
!                write(771,*) "size", ptr_Index_Intrfc(indexCommon_neighbor),&
!                     ptr_Index_Intrfc(indexCommon_neighbor+1)-1
!             End Do
!             write(771,*) "we have ", nb_of_Common_index ,"in common"
!             Do i=1,nb_of_Common_index
!                write(771,*) Array_of_Common_index(i)               
!             End Do
!          End If
!          j=1
!          Do i=ptr_Index_Intrfc(indexfailed_neighbor),&
!               ptr_Index_Intrfc(indexfailed_neighbor+1)-1
!             write(771,*)  work(index_intrfc(i))
!          End Do
         
!          write(771,*) "all x"
!          Do i=1,bound_ndof
!             write(771,*) work(i)
!          End Do
!           close(771)
!        End If
!        If (myrank .Eq. ArrayOfFaultyPorc(2)) Then
!           open(unit=772,file="Fault2.txt",action="write")   
!           write(772,*)"I am proc", myrank, "my neigbhors" 
!          Do neighb=1, nNeighb
!             write(772,*) indexVi(neighb)
!          End Do
!          if(nb_Common_neighbor > 0 ) Then
!             Do Common=1,nb_Common_neighbor
!                Do neighb=1, nNeighb
!                   If (Common_neighbor(Common) .Eq. indexVi(neighb)) Then
!                      indexCommon_Neighbor = neighb
!                      Exit  
!                   End If
!                End Do
               
!                write(772,*) "common is", Common_neighbor(Common)
!                write(772,*) "size", ptr_Index_Intrfc(indexCommon_neighbor),&
!                     ptr_Index_Intrfc(indexCommon_neighbor+1)-1
!             End Do
!             write(772,*) "we have ", nb_of_Common_index ,"in common"
!             Do i=1,nb_of_Common_index
!                write(772,*) Array_of_Common_index(i)               
!             End Do
!          End If
         
!          Do i=ptr_Index_Intrfc(indexfailed_neighbor),&
!               ptr_Index_Intrfc(indexfailed_neighbor+1)-1
!            write(772,*)  work(index_intrfc(i))
!         End Do
!         write(772,*) "all x"
!         Do i=1,bound_ndof
!            write(772,*) work(i)
!         End Do
!         close(772)
!      End If
     
9999 Continue
     If(Associated(recvbuff)) Deallocate (recvbuff)
     If(Associated(sendbuff)) Deallocate (sendbuff)
End SUBROUTINE XMPH_INTERPOLATE

SUBROUTINE MPH_COMMON_NEIGHBOR( & ! intents
       mphs,                 & ! in 
       Failedrank,           & ! in 
       Common_neighbor,      & ! out
       info                  & ! out 
       )

  Use mph_error_mod        
  Use XMPH_maphys_type
  Use MPH_maphys_enum
  Implicit None
  Include 'mpif.h'
  
  !* Arguments *!
  Type(XMPH_maphys_t)   , Intent(inout ) :: mphs
  Integer               , Intent(in  )  :: Failedrank
  Integer               , Intent(out  )  :: Common_neighbor(100)
  Integer               , Intent(out ) :: info    

    Integer :: nNeighb
    Integer :: nb_of_common_neighbor
    Integer :: neighb  
    Integer :: max_nb_neighb
    Integer :: mpiReq(2)
    Integer :: tagg 
    MPH_INT :: i,j
    ! Arrays 
    Integer  , Pointer :: indexVi       (:)
    Integer  , Pointer :: recvbuff      (:)
    Integer  , Pointer :: sendbuff      (:)
    max_nb_neighb = 100
    nb_of_common_neighbor = 0
    info =0
    tagg = 55
    nNeighb          = mphs%lc_domain%mynbvi
    indexVi          => mphs%lc_domain%myindexVi

       Allocate (recvbuff(max_nb_neighb), sendbuff(max_nb_neighb),STAT = info )
       CHCKASSRT( info == 0, info )
       If( info < 0 ) Goto 9999        
       !initialize send buff
       Do i=1,max_nb_neighb
          sendbuff(i) = 0
       End DO
       
       !fill send buf
       sendbuff(1) = nNeighb
       Do i=1,nNeighb
          sendbuff(i+1) = indexVi(i)
       End DO
       !send Data
        Call MPI_Isend(                    &
              sendbuff(1),max_nb_neighb, MPI_INT,failedrank,&
              tagg,mphs%comm, mpiReq(1),info)
         ASSRT( info == MPI_SUCCESS )       

         Call MPI_Irecv(recvbuff(1),max_nb_neighb,&
              MPI_INT,failedrank,&
              tagg,mphs%comm,mpiReq(2),info)
         ASSRT( info == MPI_SUCCESS )

         Call MPI_Waitall(2, mpiReq,MPI_STATUSES_IGNORE, info)
         ASSRT( info == MPI_SUCCESS )


         !Check if we have common neighbor
         info = 0
         Do i=1,nNeighb
            Do j=1,recvbuff(1)
               If(recvbuff(j+1) .Eq. indexVi(i)) Then
                  nb_of_common_neighbor = nb_of_common_neighbor +1
                  Common_neighbor(nb_of_common_neighbor) = indexVi(i)
               End If
            End Do
         End Do
         info = nb_of_common_neighbor
         ! set number of processors involved in recovery 
         mphs%rinfo (RINFO_PROC_INVOLVED) =recvbuff(1)+ nNeighb -info
9999 continue    
     If(Associated(recvbuff)) Deallocate (recvbuff)
     If(Associated(sendbuff)) Deallocate (sendbuff)
End SUBROUTINE MPH_COMMON_NEIGHBOR

SUBROUTINE Mph_not_contain(  & ! intents
       Array_of_Comm,         & ! in 
       nb_of_common_neighbor, & ! in 
       value,                 & ! in               
       info                   & ! out 
       )

  !* Arguments *!
  Integer          , Intent(in  )  :: Array_of_Comm(100000)
  Integer          , Intent(in )   :: nb_of_common_neighbor    
  Integer          , Intent(in )  :: value  
  Integer          , Intent(out )  :: info    

  Integer :: i

  info = 1
  If(nb_of_common_neighbor .Gt. 0) Then
      Do i=1,nb_of_common_neighbor
         If(Array_of_Comm(i) .Eq. value) Then
            info = 0
            Exit
         End If
      End Do
   End If
End SUBROUTINE Mph_not_contain
END MODULE XMPH_fault_mod

