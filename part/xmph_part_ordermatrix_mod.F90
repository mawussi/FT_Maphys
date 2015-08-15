! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"

#if ! ( HAVE_LIBSCOTCH || HAVE_LIBMETIS )     
#error "MaPHyS requires a partitioner (scotch or metis)" 
#endif

!#define MAPHYS_DEBUG 1
! [+] module : XMPH_part_order_globalmatrix_mod --------------------------
!
!> module to order the input matrix. 
!!
Module XMPH_part_ordermatrix_mod

  !* No Implicit typing *!
  Implicit None
  
  !* Access *!

  Public  :: XMPH_Part_OrderGlobalMatrix
  Private :: Part_CreateGraph
  Private :: Part_metis_partgraph
  Private :: Part_reorder_interface
  Private :: Part_tree_compute_level
  Private :: Part_tree_print
#if MAPHYS_DEBUG
  Private :: Part_CreateGraph_debug
  Private :: Part_metis_partgraph_debug
  Private :: Part_reorder_interface_debug
#endif

  !* Routines *!

  Contains

! [+] routine : XMPH_PART_OrderGlobalMatrix ------------------------------------
!
!> Part a matrix according to maphys partioning
!!
!! @param[in]     nbdom       number of domain to split into
!! @param[in]     metalgo     choosend metis algorithm
!!
!! @param[in,out] sm_global_A matrix to be splitted
!! @param[in,out] ana_timing  timings for the analysis phase
!!
!! @param[out]    graph_global_A graph associated to the matrix 
!! @param[out]    tree           partioning binary tree
!! @param[out]    info           Subroutine's status
!!
!! 
!! @par Steps :
!!
!! - [0.0] Initialize local variables
!! - [1.0] Compute the Graph 
!! - [2.0] Call METIS (Modified version)
!! - [3.0] Modify Metis' to Maphys' Permutation 
!! - [4.0] Generate the level, generate the number of the subdomains 
!! - [5.0] Save computed data, deallocate memory, etc.
!!
!! @author Yohan Lee-tin-yien
!! @todo   rename local variables
!!
Subroutine XMPH_PART_OrderGlobalMatrix    &       !intents
     ( nbdom, metalgo,        &       ! in
     sm_global_A, ana_timing, &       ! inout 
     graph_global_A, tree, info )     ! out
  
  !* Module(s) & co. *!
  Use mph_log_mod
  Use mph_error_mod
  Use MPH_part_type
  Use XMPH_sparse_matrix_type
  Use MPH_part_scotch_mod

  Implicit None

  !* Subroutine arguments *!
  Integer                     , intent(in) :: nbdom   ! number of domains 
  Integer                     , intent(in) :: metalgo ! metis algorithm

  Type(XMPH_sparse_matrix_t)  , intent(inout) :: sm_global_A
  Real(kind=8)                , intent(inout) :: ana_timing(ANA_TIMING_SIZE)

  Type(maphys_matrix_graph_t) , intent(out) :: graph_global_A
  Type(maphys_binary_tree_t)  , intent(out) :: tree
  Integer                     , intent(out) :: info

  !* External functions *!
  Real(kind=8), external :: MPI_Wtime

  !* Local Variables *!

  ! sm_global_A association
  integer          :: nnza,ndof
  MPH_INT  , pointer :: ia(:), ja(:)
  XMPH_FLOAT, pointer ::  a(:)

  ! tree association
  integer          :: toplevel
  integer, pointer :: metperm   (:)
  integer, pointer :: domstptr  (:)
  integer, pointer :: domintdof (:)
  integer, pointer :: metiperm  (:)
  integer, pointer :: domptr    (:)
  
  ! ana_timing
  real(kind=8) :: tparti
  real(kind=8) :: tmetis 
  real(kind=8) :: ttmp   

  ! other

  Integer      :: i      !< counter
  Integer      :: nproc  !< number of processus
  Integer      :: nbsep  !< number of separators
  Integer :: minprocinterior
  Integer :: maxprocinterior

  Integer      :: iinfo     !< status of calls to other Subroutines
  Integer      :: istep  
  Integer      :: msg_class !< class of log messages.
  Integer      :: verbosity !< maphys level of verbosity 
  Integer      :: logunit   !< unit used to log messages.

  Type(maphys_binary_node_t), pointer :: azz (:) !< structure used inside metis
  Character(len=MAPHYS_STRL), Parameter :: rname="PART_OrderGlobalMatrix"

  !- End of header -------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! [0] Initialize local variables
  !-----------------------------------------------------------------------------
  ! 
  nnza = sm_global_A%nnz
  ndof = sm_global_A%m
  ia   => sm_global_A%i
  ja   => sm_global_A%j
  a    => sm_global_A%v

  ! other
  nproc   = nbdom
  nbsep   = nbdom-1
  
  tparti = MPI_Wtime()
  tmetis = MPI_Wtime()

  !-----------------------------------------------------------------------------
  ! [1] Compute the Graph 
  !-----------------------------------------------------------------------------

  istep = 11
  Call Part_CreateGraph(metalgo, sm_global_A, graph_global_A , iinfo)
  If (iinfo < 0) Goto 9999

#if MAPHYS_DEBUG
  istep = 12
  Call Part_CreateGraph_debug(sm_global_A, graph_global_A, iinfo) 
  If (iinfo < 0) Goto 9999
#endif

  !-----------------------------------------------------------------------------
  ! [2] Call METIS (Modified version)
  !-----------------------------------------------------------------------------

  If ( metalgo == 4 ) Then

     !--------------------------------------------------------------------------
     ! [2.1] Call SCOTCH
     !--------------------------------------------------------------------------

     istep = 21

     Call mph_log(MSG_DEBUG, "Partitioner = SCOTCH")
     Call PART_scotch_partgraph              &  ! intents
          ( metalgo, nbdom, graph_global_A,  &  ! in
          & azz, metperm,metiperm, iinfo )       ! out
     If (iinfo < 0) Goto 9999

  Else

     !--------------------------------------------------------------------------
     ! [2.2] Call METIS
     !--------------------------------------------------------------------------
     istep = 22

     Call mph_log(MSG_DEBUG, "Partitioner = METIS")
     call Part_metis_partgraph                         & ! intents
          ( metalgo, nbdom, graph_global_A,            & ! in
          azz, metperm, metiperm, iinfo )                ! out
     If (iinfo < 0) Goto 9999

  End If
  tmetis=MPI_Wtime()-tmetis
  ttmp=MPI_Wtime()

  !--------------------------------------------------------------------------
  ! [2.3] Export to file (for verifications).
  !--------------------------------------------------------------------------

#if MAPHYS_DEBUG
  ! export permutation and matrix
  istep = 231
  Call Part_metis_partgraph_debug                    & ! intents
       (sm_global_A, metperm, metiperm, nbsep, azz, & ! in
       iinfo)                                         ! out
  If (iinfo < 0) Goto 9999
#endif

  ! print partitioning tree
  istep = 232
  msg_class = MSG_DEBUG
  Call mph_log_GetVerbosity(verbosity)
  Call mph_log_GetUnit     (msg_class, logunit)
  If ((verbosity >= msg_class).and.(logunit >= 0)) Then
     Call Part_tree_print(azz, nbsep, logunit, iinfo )
     If (iinfo < 0) Goto 9999
  End If

  !-----------------------------------------------------------------------------
  ! [3] Modify Metis' to Maphys' Permutation 
  !-----------------------------------------------------------------------------
  ! 
  ! MY NEW PERMUTATION
  ! ==================
  ! change the permutation in a way that becomes [allseparator allinterior]
  ! for that I store the new metperm in metiperm then when finished 
  ! I copy metiperm to metperm and I will generate metiperm
  !
  ! @see: [3] Maphys's Permuation (doc_ana.inc)

  istep = 31
  call Part_reorder_interface                        & 
       ( nbsep, ndof, nproc                         & ! in 
       , azz, metperm , metiperm                    & ! inout
       , domintdof , domstptr, iinfo)                 ! out     
  If (iinfo < 0) Goto 9999

  minprocinterior=minval(domintdof)
  maxprocinterior=maxval(domintdof)

  ! Permute the matrix according to MY NEW permutation
  do i=1,nnza
     ia(i)=metiperm(ia(i))
     ja(i)=metiperm(ja(i))
  enddo

#if MAPHYS_DEBUG
  istep = 32
  Call Part_reorder_interface_debug                  & ! intents
       ( sm_global_A,                               & ! in
       ndof, metperm, metiperm,                     &
       nbsep, azz,                                  &
       iinfo)                                         ! out     
  If (iinfo < 0) Goto 9999
#endif

  !-----------------------------------------------------------------------------
  ! [4] Generate the level, generate the number of the subdomains 
  !-----------------------------------------------------------------------------
  !
  ! @Warning : ONLY FOR BISSECTION GRAPH
  istep = 41
  Call Part_tree_compute_level(nbdom, nbsep, azz, domptr, toplevel, iinfo)
  If (iinfo < 0) Goto 9999

  istep = 42
  msg_class = MSG_VERBOSE
  Call mph_log_GetVerbosity(verbosity)
  Call mph_log_GetUnit     (msg_class, logunit)
  If ((verbosity >= msg_class).and.(logunit >= 0)) Then
     Call Part_tree_print(azz, nbsep, logunit, iinfo )
     If (iinfo < 0) Goto 9999
  End If

  !-----------------------------------------------------------------------------
  ! [5] Save computed data, deallocate memory, etc.
  !-----------------------------------------------------------------------------

  ! Save computed data
  tree%nbsep    =  nbsep
  tree%toplevel =  toplevel
  tree%node      => azz
  tree%domptr    => domptr
  tree%domstptr  => domstptr
  tree%domintdof => domintdof
  tree%metperm   => metperm
  tree%metiperm  => metiperm

  ana_timing(1) = tmetis
  ana_timing(2) = ttmp
  ana_timing(7) = tparti

  ! Print error/warning messages
9999 Continue
  If ( iinfo /=  0 ) Then
     
     If ( iinfo > 0) msg_class = MSG_WARNING
     If ( iinfo < 0) msg_class = MSG_ERROR
     
     Select Case(istep) 
     Case(11); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
          " while computing the graph associated to the input matrix")
     Case(12); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
          " while exporting the graph associated to the input matrix")
     Case(21); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
          " while partitioning the graph with SCOTCH")
     Case(22); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
          " while partitioning the graph with METIS")
     Case(31); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
          " while reordering the interface")
     Case(32); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
          " while checking the reordered interface")
     Case(41); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
          " while appending the level to the bissection tree")
     Case(42); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
          " while printing the bissection tree")
     End Select
     
  End If

  ! report success, failure or warning
  If ( iinfo == 0 ) info =  0
  If ( iinfo <  0 ) info = - istep
  If ( iinfo >  0 ) info = + istep

End Subroutine XMPH_PART_OrderGlobalMatrix


! [+] routine : PART_CreateGraph -----------------------------------------------
!
!> Create the graph corresponding to the matrix rows and columns 
!!
!! @param sm_global_A    the sparse matrix of the global linear system
!! @param graph_global_A the associated graph
!! @param metalgo        algorithm  based on 
!!                       - 1 = NODES, 
!!                       - 2 = EDGES,
!!                       - 3 = weighted nodes
!! @param info           status of the Subroutine (0=SUCCESS)
!!
!! @todo update algo comments
Subroutine PART_CreateGraph        &
     ( metalgo, sm_global_A,   & ! in
     & graph_global_A, info)     ! out

  !* Module(s) & co. *!
  Use XMPH_sparse_matrix_type
  Use MPH_part_type
  Use mph_log_mod
  Implicit none
         
  !* Arguments *!

  Type(XMPH_sparse_matrix_t), intent(in) :: sm_global_A
  Integer              , intent(in) :: metalgo
  Type(maphys_matrix_graph_t), intent(out) :: graph_global_A
  Integer                    , intent(out) :: info
  
  !* Local variables *!

  ! Constants
  Integer, Parameter :: UNSET= -1

  ! Scalars
  Integer               :: iinfo
  Integer               :: istep
  Integer               :: msg_class
  Integer               :: symtype
  MPH_INT :: i, j, k
  MPH_INT :: ind
  MPH_INT            :: nnza, ndof
  MPH_INT            :: maxadj
  Logical               :: TO_BE_SYMETRIZED

  ! Arrays
  MPH_INT, Pointer :: ia(:), ja(:)
  MPH_INT, Pointer :: xadj  (:)   ! xadjacency
  MPH_INT, Pointer :: adjncy(:)   ! adjacency 
  MPH_INT, Pointer :: vwgt  (:)   ! weight on vertices

  ! Strings
  Character(len=MAPHYS_STRL) :: rname="PART_CreateGraph"

  !- End of header ------------------------------------------------------------

  !----------------------------------------------------------------------------
  ! [1] Initialize and associate local variables
  !----------------------------------------------------------------------------

  ! init graph
  graph_global_A%ndof   = UNSET
  graph_global_A%nnz    = UNSET
  graph_global_A%maxadj = UNSET
  graph_global_A%algo   = UNSET

  Nullify(graph_global_A%xadj)
  Nullify(graph_global_A%adjncy)
  Nullify(graph_global_A%vwgt)

  ! init local scalars
  
  iinfo = 0
  symtype = sm_global_A%sym
  nnza = sm_global_A%nnz
  ndof = sm_global_A%m
  ia   => sm_global_A%i
  ja   => sm_global_A%j
  

  Select Case(symtype)
  Case(SM_SYM_IsGeneral  ) 
     TO_BE_SYMETRIZED = .False. ! structure was already symmetrized
  Case(SM_SYM_IsSymmetric,SM_SYM_IsSPD) 
     TO_BE_SYMETRIZED = .True.  
  Case Default            
     TO_BE_SYMETRIZED = .True.  
  End Select

  !----------------------------------------------------------------------------
  ! [2] Compute the graph according to its symmetry
  !----------------------------------------------------------------------------
  If( TO_BE_SYMETRIZED )Then

     !-------------------------------------------------------------------------
     ! [2.1] Transform (ia,ja) into graph data structure
     ! [ - ] (for symmetric --> triangular entries)
     !------------------------------------------------------------------------- 

     !-------------------------------------------------------------------------
     ! [2.1.1] Allocate xadj,adjcncy, vwgt
     !------------------------------------------------------------------------- 

     istep = 211
     Allocate(xadj(ndof+1),STAT=iinfo)
     IF( iinfo /= 0) iinfo = -iinfo
     If (iinfo /= 0) Goto 9999

     Allocate(adjncy(2*(nnza)),STAT=iinfo)
     IF( iinfo /= 0) iinfo = -iinfo
     If (iinfo /= 0) Goto 9999

     If(metalgo == 3) Then
        Allocate(vwgt(ndof),STAT=iinfo)
        IF( iinfo /= 0) iinfo = -iinfo
        If (iinfo /= 0) Goto 9999
     End If

     !-------------------------------------------------------------------------
     ! [2.1.2] Construct xadj 
     !------------------------------------------------------------------------- 

     ! init
     Do i=1,ndof+1
        xadj(i) = 0
     End Do

     ! count the number of edges associated to each row 
     xadj(1)=1
     do k=1,nnza
        i=ia(k)
        j=ja(k)
        if(i /= j)then
           xadj(i+1)= xadj(i+1) + 1
           xadj(j+1)= xadj(j+1) + 1
        endif
     enddo

     ! sum the counts to form xadj
     xadj(1)=1
     do i=2,ndof+1
        xadj(i)= xadj(i) + xadj(i-1)
     enddo

     !-------------------------------------------------------------------------
     ! [2.1.3] Form the adjncy
     !------------------------------------------------------------------------- 
     
     !> @warning - about xadj
     !! xadj is complete at this step.
     !! but momentarily modify it to save memory. 
     !! value is retrieve after shifting.

     ! construct adjcncy (& modify xadj)
     Do k=1,nnza
        i=ia(k)
        j=ja(k)
        If(i /= j)Then

           ind          = xadj(i)
           adjncy(ind)  = j
           xadj(i)      = xadj(i) + 1 

           ind          = xadj(j)
           adjncy(ind)  = i
           xadj(j)      = xadj(j) + 1 

        End If
     End Do

     ! retrieve xadj by shifting it by 1
     Do k=1,ndof
        xadj(ndof-k+2) = xadj(ndof-k+1)
     End Do
     xadj(1)=1

     !-------------------------------------------------------------------------
     ! [2.1.4] Form the weight of vertices (here weight=adjncy)  
     !------------------------------------------------------------------------- 

     If (metalgo == 3) Then
        Do i=1,ndof
           vwgt(i) = xadj(i+1)-xadj(i)
        End Do
     Endif

  Else ! NOT TO_BE_SYMMETRIZED

     !-------------------------------------------------------------------------
     ! [2.2] transform (ia,ja) into graph data structure
     ! [ - ] (for symmetric pattern matrices --> all entries)
     !-------------------------------------------------------------------------
     ! See also section [2.1] for full comments

     ! allocate xadj, adjcny and vwgt 
     istep = 221
     Allocate(xadj(ndof+1),STAT=iinfo)
     IF( iinfo /= 0) iinfo = -iinfo
     If (iinfo /= 0) Goto 9999

     Allocate(adjncy(nnza),STAT=iinfo)
     IF( iinfo /= 0) iinfo = -iinfo
     If (iinfo /= 0) Goto 9999

     If(metalgo == 3) Then
        Allocate(vwgt(ndof),STAT=iinfo)
        IF( iinfo /= 0) iinfo = -iinfo
        If (iinfo /= 0) Goto 9999
     End If

     ! form xadj
     Do i=1,ndof+1
        xadj(i) = 0
     End Do

     xadj(1)=1
     do k=1,nnza
        i=ia(k)
        j=ja(k)
        if(i /= j)then
           xadj(i+1)= xadj(i+1) + 1
        endif
     enddo

     xadj(1)=1
     do i=2,ndof+1
        xadj(i)= xadj(i) + xadj(i-1)
     enddo

     ! form the adjncy
     do k=1,nnza
        i=ia(k)
        j=ja(k)
        if(i /= j)then
           ind          = xadj(i)
           adjncy(ind)  = j
           xadj(i)    = xadj(i) + 1 
        endif
     enddo

     do i=1,ndof
        xadj(ndof-i+2) = xadj(ndof-i+1)
     enddo
     xadj(1)=1

     ! form the weight of vertices (here weight=adjncy)  
     If(metalgo == 3) Then
        do i=1,ndof
           vwgt(i) = xadj(i+1)-xadj(i)
        enddo
     Endif

  end IF

  !-----------------------------------------------------------------------------
  ! [3] Compute statistics
  !-----------------------------------------------------------------------------

  maxadj = 0
  Do i=1,ndof
     k= xadj(i+1)-xadj(i) 
     maxadj= Max(maxadj,k)
  End Do
  
  !-----------------------------------------------------------------------------
  ! [4] Finish
  !-----------------------------------------------------------------------------

  ! Save computed data
  graph_global_A%ndof  = ndof
  graph_global_A%nnz   = nnza
  graph_global_A%maxadj   = maxadj

  graph_global_A%xadj    => xadj
  graph_global_A%adjncy  => adjncy

  graph_global_A%algo  = metalgo
  If ( metalgo == 3) graph_global_A%vwgt => vwgt

  ! Print error messages
9999 Continue
  If ( iinfo /=  0 ) Then

     If ( iinfo > 0) msg_class = MSG_WARNING
     If ( iinfo < 0) msg_class = MSG_ERROR
     
     Select Case(istep) 
     Case(211); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
          " while allocating data for the graph")
     Case(221); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
          " while allocating data for the graph")
     End Select
     
  End If

  ! report success, failure
  If ( iinfo == 0 ) info =  0
  If ( iinfo <  0 ) info = - istep
  If ( iinfo >  0 ) info =   istep

End Subroutine Part_CreateGraph

#if MAPHYS_DEBUG
! [+] routine : PART_CreateGraph_debug -----------------------------------------
!
!> Check routine PART_CreateGraph
!!
!! Export the graph and the matrix to files if MAPHYS_DEBUG is defined and positive.
!! If not, do nothing.
!! 
!! @param[in]  sm_global_A     The input matrix to write to file
!! @param[in]  graph_global_A  The graph to write to file
!! @param[out] info            The routine status 
!!
!! @author Azzam Haidar
!! @author Yohan Lee-tin-yien
!!
Subroutine PART_CreateGraph_debug        &
     ( sm_global_A, graph_global_A,  & ! in
     & info)                           ! out

  !* Modules *!
  Use MPH_part_type
  Use XMPH_sparse_matrix_mod
  Use mph_log_mod
  Use mph_dbg_mod
  Implicit None

  !* Arguments *!
  Type(XMPH_sparse_matrix_t)      , intent(in)  :: sm_global_A
  Type(maphys_matrix_graph_t), intent(in)  :: graph_global_A
  Integer                    , intent(out) :: info
  
  !* local variables *!
  ! Scalars
  Integer :: i, j, ndof
  Integer :: iinfo

  ! Arrays
  MPH_INT, pointer :: xadj  (:)   ! xadgency
  MPH_INT, pointer :: adjncy(:)   ! adgency 

  ! End of header --------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! [1] Export the graph
  !-----------------------------------------------------------------------------

  Call MPH_dbg_init

  Call MPH_dbg_set_file("gb_A_Graph.txt")
  xadj    => graph_global_A%xadj
  adjncy  => graph_global_A%adjncy
  ndof = sm_global_A%m
  Do i=1,ndof
     Write(dbg_unit,*)"xadj  :",i,xadj(i),xadj(i+1)-1,"   : ",&
          (adjncy(j),j=xadj(i),xadj(i+1)-1)
     Write(dbg_unit,*)
  Enddo

  !-----------------------------------------------------------------------------
  ! [2] Export the matrix
  !-----------------------------------------------------------------------------

  Call MPH_dbg_set_file("Part_CreateGraph_gb_A.mtx")
  Call XMPH_sm_mmwrite(sm_global_A,dbg_unit, iinfo) 

  !-----------------------------------------------------------------------------
  ! [3] Report status
  !-----------------------------------------------------------------------------

  info = iinfo

end Subroutine Part_CreateGraph_debug
#endif  


! [+] routine : Part_metis_partgraph --------------------------------------------
!
!> partition the matrix graph with metis
!!
!! @param [in ] metalgo     choosen algorithm
!! @param [in ] nbdom       wanted number of domain
!! @param [out] azz         the binary tree (size = nbdom -1)
!! @param [out] metperm     the permutation
!! @param [out] metiperm    the inverse permutation
!! @param [out] info        status of the Subroutine (0=SUCCESS)
!!
Subroutine Part_metis_partgraph              &  ! intents
     ( metalgo, nbdom, graph_global_A,      &  ! in
     & azz, metperm,metiperm, info )           ! out

  !* Module(s) *!
  Use MPH_part_type
  Use mph_log_mod
  Implicit None

  !* Arguments *!
  Integer   , intent(in) :: metalgo
  Integer   , intent(in) :: nbdom
  Type(maphys_matrix_graph_t), intent(in) :: graph_global_A
  Type(maphys_binary_node_t) , pointer, intent(  out) :: azz (:)
  Integer, intent(out), pointer :: metperm  (:)
  Integer, intent(out), pointer :: metiperm (:)
  Integer, intent(out)          :: info

#ifndef HAVE_LIBMETIS
  Call MPH_Log(MSG_ERROR,"MaPHyS was not compiled with METIS (-DHAVE_LIBMETIS)")
  info = -1
  Return
#else
  !* local variables *!

  ! constants
  Integer , Parameter :: noptions = 11
  integer, parameter :: BEGIN=1, SUCCESS=0, FAIL= -1 

  ! scalars
  Integer :: ndof
  Integer :: i
  Integer :: numflag
  Integer :: iinfo, istep, msg_class

  ! arrays
  integer  :: options(noptions)
  MPH_INT, pointer :: xadj(:) , adjncy (:) , vwgt(:)

  ! String 
  Character(len=MAPHYS_STRL) :: rname = "PART_METIS_PARTGRAPH"

  !- End of header -------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! [0] Initialize local data
  !-----------------------------------------------------------------------------

  iinfo = SUCCESS

  ndof    =  graph_global_A%ndof
  xadj    => graph_global_A%xadj
  adjncy  => graph_global_A%adjncy
  if ( metalgo == 3) vwgt => graph_global_A%vwgt
  
  istep = 01
  If( nbdom < 1 ) iinfo = FAIL
  If( iinfo /= SUCCESS ) Goto 9999

  istep = 02
  ALLOCATE(azz(nbdom-1), stat=iinfo )
  IF( iinfo /= SUCCESS ) iinfo = -iinfo
  IF( iinfo /= SUCCESS ) Goto 9999

  Do i=1,nbdom-1
     azz(i)%lged   = 0
     azz(i)%father = 0
     azz(i)%rson   = 0
     azz(i)%lson   = 0
     azz(i)%sepst  = 0
     azz(i)%seped  = 0
     azz(i)%rgst   = 0
     azz(i)%rged   = 0
     azz(i)%lgst   = 0
     azz(i)%level  = 0
     azz(i)%rdomid = -1
     azz(i)%ldomid = -1
  End Do

  !-----------------------------------------------------------------------------
  ! [1] Call metis
  !-----------------------------------------------------------------------------

  numflag = 1
  Do i=1, noptions
     options(i)= 0
  End Do

  istep = 11
  ALLOCATE(metperm(ndof), stat=iinfo )
  IF( iinfo /= SUCCESS ) iinfo = -iinfo 
  IF( iinfo /= SUCCESS ) Goto 9999

  istep = 12
  ALLOCATE(metiperm(ndof), stat=iinfo )
  IF( iinfo /= SUCCESS ) iinfo = -iinfo
  IF( iinfo /= SUCCESS ) Goto 9999


  SELECT CASE(metalgo)
  CASE(1)
     options(1)= 1   ! (0 use default 1 use MY options)
     options(2)= 3   ! OPTION_CTYPE = MATCH_SHEM=3
     options(3)= 1   ! OPTION_ITYPE = IPART_GGPKL=1
     options(4)= 2   ! OPTION_RTYPE = RTYPE_SEP1SIDED=2

     options(5)= 0   ! OPTION_DBGLVL= 0
     options(6)= 1   ! OPTION_OFLAGS= OFLAG_COMPRESS=1  (compress flag 0 or 1 (default))
     options(7)= -1  ! OPTION_PFACTOR = -1
     options(8)= 5   ! OPTION_NSEPS = 1

     CALL PO_METIS_NodeND (ndof,xadj,adjncy,numflag, &
     &               options,metperm,metiperm,nbdom,azz)
  CASE(2)
     options(1)= 0 
     CALL PO_METIS_EdgeND (ndof,xadj,adjncy,numflag, &
     &               options,metperm,metiperm,nbdom,azz)

  CASE(3)
     options(1)= 0 

     CALL PO_METIS_NodeWND(ndof,xadj,adjncy,vwgt,numflag, &
     &               options,metperm,metiperm,nbdom,azz)
  CASE DEFAULT
     istep = 131
     iinfo = FAIL
     If (iinfo /= SUCCESS) Goto 9999
  END SELECT

  !-----------------------------------------------------------------------------
  ! [2] Handle Errors
  !-----------------------------------------------------------------------------

  ! Print error/warning messages
9999 Continue
  If ( iinfo /=  0 ) Then
     
     If ( iinfo > 0) msg_class = MSG_WARNING
     If ( iinfo < 0) msg_class = MSG_ERROR
     
     Select Case(istep) 
     Case(01); Call mph_logWithInfo (msg_class,nbdom,Trim(rname)//&
          "Bad value of nbdom =")
     Case(02); Call mph_logWithInfo (msg_class,nbdom-1,Trim(rname)//&
          "Failing to allocate the leafs of the paritioning tree, size=")
     Case(11 ); Call mph_logWithInfo (msg_class,ndof,Trim(rname)//&
          "Failing to allocate metperm, size=")
     Case(12 ); Call mph_logWithInfo (msg_class,ndof,Trim(rname)//&
          "Failing to allocate metiperm, size=")
     Case(131); Call mph_logWithInfo (msg_class,metalgo,Trim(rname)//&
          "Bad partitioning strategy for METIS, value=")
     End Select
     
  End If

  ! report success, failure or warning
  If ( iinfo == 0 ) info =  0
  If ( iinfo <  0 ) info = - istep
  If ( iinfo >  0 ) info = + istep
#endif

end Subroutine Part_metis_partgraph

!
!-------------------------------------------------------------------------------
!
#if MAPHYS_DEBUG
Subroutine Part_metis_partgraph_debug&
     (sm_global_A, metperm, metiperm,&
     nbsep, azz, info)
  
  !* Module(s) *!
  Use MPH_part_type
  Use XMPH_sparse_matrix_mod
  Use mph_dbg_mod
  Implicit none

  !* Arguments *!

  type(XMPH_sparse_matrix_t)       , intent(in) :: sm_global_A
  integer, intent(in), pointer            :: metperm  (:)
  integer, intent(in), pointer            :: metiperm (:)

  integer                     , intent(in) :: nbsep
  type(maphys_binary_node_t), pointer, intent(in) :: azz (:)
  integer              , intent(out)  :: info

  !* Local variables *!

  MPH_INT :: i
  MPH_INT :: n
  integer :: iinfo

  !- End of header -------------------------------------------------------------

  !-----------------------------------------------------------------------------   
  ! [1] Init
  !-----------------------------------------------------------------------------   

  n = sm_global_A%n

  !-----------------------------------------------------------------------------   
  ! [2] Write permutations
  !-----------------------------------------------------------------------------   

  Call MPH_dbg_set_file("metperm.txt")
  Do i=1,n
     Write(dbg_unit,*) metperm(i)
  End Do

  Call MPH_dbg_set_file("metiperm.txt")
  Do i=1,n
     Write(dbg_unit,*) metperm(i)
  End Do

  !-----------------------------------------------------------------------------   
  ! [3] Write matrix
  !-----------------------------------------------------------------------------   

  Call MPH_dbg_set_file("gb_A_afterPartionerCall.mtx")
  Call sparse_matrix_mmwrite(sm_global_A, dbg_unit, iinfo )

  !-----------------------------------------------------------------------------   
  ! [4] Write the bissection tree
  !-----------------------------------------------------------------------------   
  
  Call MPH_dbg_set_file("ptree.txt")
  Call Part_tree_print(azz, nbsep, dbg_unit, iinfo )

  info = 0

end Subroutine Part_metis_partgraph_debug

#endif



!
! [+] routine : Part_reorder_interface ------------------------------------------
!
!> modify the metis partition (== permutations) for maphys
!!
!! on output : metperm and perm generate the following ordering
!! Matrix is ordered as follow
!!         col                    col
!! row [separator_row       connex_sep-int]
!! row [interior__row       connex_int-sep]
!! that is:        
!! separator_row [ Abb  Abi1   Abi2    Abi3   Abi4   Abi5]
!! interior__row [ Aib1 Aii1    0       0      0       0 ]
!! interior__row [ Aib2   0    Aii2     0      0       0 ]
!! interior__row [ Aib3   0     0      Aii3    0       0 ]
!! interior__row [ Aib4   0     0       0     Aii4     0 ]
!! interior__row [ Aib5   0     0       0      0     Aii5]
!!
Subroutine Part_reorder_interface &
     ( nbsep, ndof, nproc        & 
     , azz, metperm , metiperm   &
     , domintdof , domstptr, info)
  Use MPH_part_type
  Implicit None
  ! IN
  integer :: nbsep
  integer :: ndof
  integer :: nproc

  ! INOUT
  MPH_INT  , pointer, intent(inout) :: metperm  (:)
  MPH_INT  , pointer, intent(inout) :: metiperm (:)
  type(maphys_binary_node_t)  , intent(inout) :: azz  (:)

  ! OUT 
  integer, pointer, intent(out) :: domintdof (:)
  integer, pointer, intent(out) :: domstptr  (:)
  integer, intent(out) :: info

  ! local variables
  integer :: totinterface
  integer :: i, j,  k
  integer :: domnb
  integer :: domst
  integer :: sizdom
  integer :: separatorst

  !     ==============================================================
  !                       MY NEW PERMUTATION
  !     ==============================================================
  !     change the permutation in a way that becomes [allseparator allinterior]
  !     for that I store the new metperm in metiperm then when finished 
  !     I copy metiperm to metperm and I will generate metiperm
  !     ==============================================================

  ALLOCATE(domintdof(nproc),STAT=info)
  If (info /= 0) info = -1
  If (info /= 0) Return

  ALLOCATE(domstptr(nproc),STAT=info)
  If (info /= 0) info = -1
  If (info /= 0) Return

  Do i=1,nproc
     domintdof(i)=0
     domstptr(i)=0
  End Do

  totinterface = 0
  do i=1, nbsep 
     totinterface=totinterface+(azz(i)%seped-azz(i)%sepst+1)
  enddo

  separatorst = 1
  domst       = totinterface+1
  domnb=-1
  do i=1, nbsep 
     !----------------------
     ! interior permutation
     !----------------------
     if(azz(i)%rson.eq.-1) then 
        k=domst-1
        do j=azz(i)%rgst,azz(i)%rged
           k=k+1
           metiperm(k)=metperm(j)
        enddo
        sizdom=azz(i)%rged-azz(i)%rgst+1
        domnb=domnb+1
        domintdof(domnb+1)= sizdom
        domstptr (domnb+1)= domst
        azz(i)%rgst = domst
        azz(i)%rged = domst+sizdom-1
        domst = domst+sizdom
        k=domst-1
        do j=azz(i)%lgst,azz(i)%lged
           k=k+1
           metiperm(k)=metperm(j)
        enddo
        sizdom=azz(i)%lged-azz(i)%lgst+1
        domnb=domnb+1
        domintdof(domnb+1)= sizdom
        domstptr (domnb+1)= domst
        azz(i)%lgst = domst
        azz(i)%lged = domst+sizdom-1
        domst = domst+sizdom
     else !
        azz(i)%rgst = 0
        azz(i)%rged = 0
        azz(i)%lgst = 0
        azz(i)%lged = 0
     endif
     !----------------------
     ! separator permutation
     !----------------------
     k=separatorst-1
     do j=azz(i)%sepst,azz(i)%seped
        k=k+1
        metiperm(k)=metperm(j)
     enddo
     sizdom=azz(i)%seped-azz(i)%sepst+1
     azz(i)%sepst = separatorst
     azz(i)%seped = separatorst+sizdom-1
     separatorst  = separatorst+sizdom
  enddo

  ! copy metiperm to metperm
  Do i=1, ndof
     metperm(i) = metiperm(i)
  End Do

  ! create the iperm
  do i=1,ndof
     metiperm( metperm(i) )=i
  enddo

end Subroutine Part_reorder_interface

!
!-------------------------------------------------------------------------------
!
#if MAPHYS_DEBUG
Subroutine Part_reorder_interface_debug & ! intents
     ( sm_global_A,                    & ! in
     & ndof, metperm, metiperm,        &
     & nbsep, azz,                     &
     & info)                             ! out     

  !* Modules *!
  Use mph_dbg_mod
  Use XMPH_sparse_matrix_mod
  Use MPH_part_type
  implicit none

  !* Arguments *!

  ! in
  type(XMPH_sparse_matrix_t), intent(in ) :: sm_global_A

  integer, intent(in)          :: ndof
  integer, intent(in), pointer :: metperm  (:)
  integer, intent(in), pointer :: metiperm (:)

  integer, intent(in) :: nbsep
  type(maphys_binary_node_t), pointer, intent(in) :: azz (:)

  ! out
  integer              , intent(out) :: info

  !* Local variables *!

  ! constants
  integer, parameter :: SUCCESS=0, ERROR=-1 

  ! scalars
  integer :: ios1
  MPH_INT :: i, j, k

  ! types
  type(XMPH_sparse_matrix_t) :: sm_ori

  ! End of header --------------------------------------------------------------



  Call MPH_dbg_init

  Call MPH_dbg_set_file("gb_A_permutated.mtx")
  Call sparse_matrix_mmwrite(sm_global_A, dbg_unit, ios1)


  Call MPH_dbg_set_file("metperm.txt")
  Do i=1, ndof 
     Write(unit=dbg_unit,FMT='(I16)') metperm(i)
  Enddo
     

  Call MPH_dbg_set_file("gb_A_reconstructed.mtx")
  Call sparse_matrix_dup(sm_ori,sm_global_A,ios1)
  Do k = 1, sm_ori%nnz
     sm_ori%i(k) = metperm(sm_global_A%i(k))
     sm_ori%j(k) = metperm(sm_global_A%j(k))
  End Do
  Call sparse_matrix_mmwrite(sm_ori, dbg_unit, ios1)


  Call MPH_dbg_set_file("separator.txt")
  Do i=1, nbsep 
     Do j=azz(i)%sepst,azz(i)%seped
        Write(UNIT=dbg_unit,FMT='(I16)')metperm(j)
     Enddo
  Enddo
  
  Call MPH_dbg_set_file("metiperm.txt")
  Do i=1, ndof 
     Write(UNIT=dbg_unit,FMT='(I16)') metiperm(i) 
  End Do

  info = SUCCESS

End Subroutine  Part_reorder_interface_debug
#endif


! [+] routine : part_tree_compute_level ----------------------------------------
!
!> Append the level information to a partitioning binary tree
!!
!! Append to each separator in the partitioning binary tree
!! its level in the binary tree.
!! Also compute the level of the root, "toplevel" (= the tree depth).
!! And "domptr", such as : "domptr(domidx) = sepidx"
!! with domidx the index of a domain    ( 1..nbdom )
!! and sepidx  the index of a separator ( 1..nbsep )
!!
!!-----
!! 
!! @param [in]      nbdom     number of domain in the partition
!! @param [in]      nbsep     number of separators int the tree
!! @param [in,out]  azz       list of separators
!! @param [   out]  domptr    see above paragraph
!! @param [   out]  toplevel  depth of the tree
!! @param [   out]  info      routine status
!!
!!-----
!!
!! @author Yohan Lee-tin-yien
!! @author Azzam Haidar
!!
Subroutine PART_tree_compute_level(nbdom, nbsep, azz, domptr, toplevel, info)

  !* Module(s) *!
  Use MPH_part_type
  Use mph_log_mod
  Implicit None

  !* Arguments *!
  Integer, intent(in) :: nbdom
  Integer, intent(in) :: nbsep
  Type(maphys_binary_node_t), pointer, intent(inout) :: azz (:)
  Integer, pointer, intent(out) :: domptr   (:)
  Integer         , intent(out) :: toplevel
  Integer         , intent(out) :: info

  !* Local variables *!

  ! Scalars
  Integer :: i, j, k
  Integer :: iinfo,istep, msg_class
  Integer :: levid, domnb, fatherid, sepid, newlev

  ! Strings
  Character(len=MAPHYS_STRL), Parameter :: rname = "PART_tree_compute_level"
  Character(len=MAPHYS_STRL) :: msg

  !- End of header -------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! [1] Init
  !-----------------------------------------------------------------------------

  ! allocate memory
  ! domptr is a pointer that point the structure sepid of each domain
  istep = 11
  ALLOCATE(domptr(nbdom),STAT=iinfo)
  If (iinfo /= 0 ) iinfo = -iinfo
  If (iinfo /= 0 ) Goto 9999

  ! force computation of levels
  Do i= 1, nbsep
     azz(i)%level=0
  End Do

  ! init counters
  k=0
  levid = 0
  domnb = -1 ! 0 pour debuter a 1 ou -1 pour debuter a zero confondu avec procid 
  toplevel = 0 
  levid=levid+1

  !-----------------------------------------------------------------------------
  ! [2] Compute 
  !-----------------------------------------------------------------------------

  istep = 2
  do i=1, nbsep 
     if(azz(i)%rson.eq.-1) then 
        azz(i)%level  = levid

        ! right
        k=k+1
        domptr(k)=i
        domnb         = domnb+1
        azz(i)%rdomid = domnb

        ! left
        k=k+1
        domptr(k)=i
        domnb         = domnb+1
        azz(i)%ldomid = domnb

        sepid         = i
        newlev        = levid
        do j=1,nbdom-1

           ! test the end of tree
           fatherid = azz(sepid)%father
           if(fatherid.eq.0)  then 
              toplevel = azz(sepid)%level
              goto 319 ! end of tree go out
           endif

           ! append the level
           newlev=newlev+1
           if(azz(fatherid)%level.eq.0) then ! do not yet have a level 
              azz(fatherid)%level  = newlev
           else 

              if (azz(fatherid)%level /= newlev)then
                 Write(msg,'(2A,I5,I5)') Trim(rname)," different number of level found",&
                      azz(fatherid)%level,newlev
                 Call mph_log(MSG_ERROR,msg)
                 iinfo = -1
                 Goto 9999
              endif

              goto 319
           endif
           sepid    = fatherid
        enddo
     endif
319  continue
  enddo

  !-----------------------------------------------------------------------------
  ! [3] Finish
  !-----------------------------------------------------------------------------

  ! print messages
9999 Continue
  If ( iinfo /=  0 ) Then
     
     If ( iinfo > 0) msg_class = MSG_WARNING
     If ( iinfo < 0) msg_class = MSG_ERROR
     
     Select Case(istep) 
     Case(11); Call mph_logWithInfo (msg_class,nbdom,Trim(rname)//&
          " failed to allocate memory, domptr, size=")
     Case( 2); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
          " internal error in the bissection tree")
     End Select
     
  End If
  
  ! report success, failure or warning
  If ( iinfo == 0 ) info =  0
  If ( iinfo <  0 ) info = -istep
  If ( iinfo >  0 ) info =  istep

end Subroutine Part_tree_compute_level


! [+] routine : part_tree_print -----------------------------------------------
!
!> Print the partitioning binary tree.
!!
!! @param [in ] azz    The leafs of the binary tree
!! @param [in ] nbsep  The number of separator in the binary tree
!! @param [in ] wunit   The unit to write the data.
!! @param [out] info   The routine status
!!
!! @author Yohan Lee-tin-yien
!! @author Azzam Haidar
!!
Subroutine Part_tree_print (azz, nbsep, wunit, info)

  !* Module(s) & co. *!
  Use MPH_part_type
  Implicit none

  !* Arguments *!

  Type(maphys_binary_node_t), pointer, intent(in) :: azz (:)
  Integer , intent(in)                            :: nbsep
  Integer , intent(in)                            :: wunit
  Integer ,intent(out)                            :: info 

  !* Local variables *!
  Integer :: i
  Character(len=MAPHYS_STRL) :: gra_print

  !- End of header -------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! [1] Init
  !-----------------------------------------------------------------------------

  Write(gra_print,'(A)') '==========================================&
       &==================================================================&
       &=========================' 
  
  !-----------------------------------------------------------------------------
  ! [2] Write the tabular
  !-----------------------------------------------------------------------------
  
  ! header
  Write(wunit,'(A)') Trim(gra_print)
  Write(wunit,'(50X,A)') '= Partitioning bissection tree ='
  Write(wunit,'(13A10)') "Id","father","rson","lson","sepst","seped",   &
       "rgst","rged","lgst","lged",                       &
       'rdomid','ldomid','level'
  Write(wunit,'(A)') Trim(gra_print)
  
  ! data
  Do i=1,nbsep
     Write(wunit,'(13I10)') i, azz(i)%father, azz(i)%rson, azz(i)%lson, &
          azz(i)%sepst, azz(i)%seped, azz(i)%rgst,            &
          azz(i)%rged, azz(i)%lgst, azz(i)%lged,              &
          azz(i)%rdomid, azz(i)%ldomid, azz(i)%level
  Enddo
  
  !
  Write(wunit,'(A)') Trim(gra_print)
  
  !-----------------------------------------------------------------------------
  ! [3] Finish
  !-----------------------------------------------------------------------------

  info = 0

End Subroutine Part_tree_print



End Module XMPH_part_ordermatrix_mod
