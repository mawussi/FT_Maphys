! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"

! [+] module : XMPH_part_builddomains_mod --------------------------------------
!
!> module to build the domains (ie the interfaces and the interiors).
!!
Module XMPH_part_builddomains_mod

  !* Modules *!
  Use MPH_error_mod

  !* No Implicit typing *!
  Implicit None
  
  !* Private constants *!
  Character(len=MAPHYS_STRL), Parameter, Private :: &
       FLNAME = "XMPH_ARITHpart_builddomains_mod.F90"
  
  !* Access *!

  Public  ::  XMPH_Part_build_domains
  Private ::  XMPH_Part_split_interface
  Private ::  combivtxpcsz_realloc
  Private ::  vtxpcsz_realloc
  Private ::  XMPH_Part_sort_GlobalMatrix_CSR 
  Private ::  XMPH_Part_compute_lagrange 
#if MAPHYS_DEBUG
  Private ::  XMPH_Part_compute_lagrange_debug 
#endif
  Private ::  XMPH_Part_compute_domains_nnz 

  !* Routines *!

  Contains

! [+] routine : XMPH_part_builddomains_mod ---------------------------------------
!
!> Compute all domains (essentially the interfaces) 
!! Identify all domains and performs some treatement on them.
!!
!! @author Yohan Lee-tin-yien
!! 
!! @par Steps :
!! - [0.0] Initialize local variables
!! - [1.0] Split the interface
!! - [2.0] Sort the matrix into CSR + detect rows/columns that needs Lagrange
!! - [3.0] Apply Lagrange treatement
!! - [4.0] Compute nnzs for each domain, interface, interior
!! - [5.0] Save results, free memory, etc.
!!
!!----
!!
!! @param [in,out] sm_global_A
!!
!!      The input matrix.
!!      - On input, it must be previously sorted such as :
!!        - the first bloc of rows corresponds to the interface
!!        - the other blocs "i" of rows corresponds 
!!          to the interiors of domain "i-1"
!!        - this special sort was done in part_order_global_matrix().
!!
!!      - On output, the matrix is sorted and in CSR format.
!!
!! @param [in,out] graph
!!
!!      The graph associated to the input matrix.
!!
!! @param [in,out] tree
!!
!!      The nested dissection binary tree associated. 
!!      - On output, the structure is freed of all its contents.
!!        Several of its fields are passed to the structure all_domains.
!!
!! @param [in,out] ana_timing
!!
!!      An array holding the durations of treatements in the analysis step.  
!!      - On input, it contains the previous timings
!!      - On output, we append the durations :
!!        - to treat the data from the partitioner and 
!!          form MaPHyS Domain Decomposition.
!!        - to sort the matrix into CSR
!!
!! @param [   out] all_domains
!!
!!      Structure specifiing the MaPHyS domain decomposition.
!!     
!! @param [   out] info
!! 
!!      The routine status
!!
!!----
!! 
!! @author Azzam Haidar
!! @author Yohan Lee-tin-yien
!!
Subroutine XMPH_Part_build_domains &
     (sm_global_A, graph, tree, ana_timing, &
     all_domains, info )
  
  !* Module(s) *!

  Use MPH_part_type
  Use XMPH_sparse_matrix_mod
  Use mph_log_mod
  Implicit None

  !* Arguments *!
  
  Type(XMPH_sparse_matrix_t) , intent(inout) :: sm_global_A
  Type(maphys_matrix_graph_t), intent(inout) :: graph 
  Type(maphys_binary_tree_t) , intent(inout) :: tree
  Real(kind=8), intent(inout) :: ana_timing(ANA_TIMING_SIZE)

  Type(maphys_domains_t) , intent(out) :: all_domains
  Integer                , intent(out) :: info  

  !* External functions *!
  
  Real(kind=8), external :: MPI_Wtime

  !* Local Variables *!

  ! Constants 
  Character(len=MAPHYS_STRL), Parameter :: rname="PART_Update_partitioning"

  ! Scalars
  Integer           :: iinfo, istep, msg_class ! error/log related
  MPH_INT        :: nnza ! matrix number of entries
  Integer           :: ndof ! matrix order
  Real(kind=8)      :: ttmp  ! timing 
  Real(kind=8)      :: tsort ! timing

  ! Arrays
  Integer, Pointer  :: domptr   (:) ! partitioning tree related, pointer to domains
  Integer, Pointer  :: metperm  (:) ! partitioning tree related, permutation
  Integer, Pointer  :: domLg    (:) ! lagrange related
  Integer, Pointer  :: need_Lg (:)  ! lagrange related

  !- End of header -------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! [0] Initialize local variables
  !-----------------------------------------------------------------------------

  ndof = sm_global_A%m
  nnza = sm_global_A%nnz
  
  domptr  => tree%domptr
  metperm => tree%metperm

  ttmp    = ana_timing(2)

  !-----------------------------------------------------------------------------
  ! [1] Split the interfaces
  !-----------------------------------------------------------------------------

  istep = 1
  Call  XMPH_Part_split_interface(graph, tree, all_domains,iinfo)
  If (iinfo < 0 ) Goto 9999

  !-----------------------------------------------------------------------------
  ! [2] Sort the matrix into CSR + detect rows/columns that needs Lagrange
  !-----------------------------------------------------------------------------

  istep = 21
  Allocate(need_Lg(nnza),STAT=iinfo)
  If (iinfo >  0 ) iinfo = -iinfo
  If (iinfo <  0 ) Goto 9999

  istep = 22
  Call XMPH_Part_sort_GlobalMatrix_CSR( sm_global_A, need_Lg, tsort, iinfo )
  If (iinfo < 0 ) Goto 9999
  
  !-----------------------------------------------------------------------------
  ! [3] Apply Lagrange treatement
  !-----------------------------------------------------------------------------
  
  istep = 31
  Call XMPH_Part_compute_lagrange (tree, need_Lg, domLg, iinfo )
  If (iinfo < 0 ) Goto 9999

#if MAPHYS_DEBUG
  istep = 32
  Call XMPH_Part_compute_lagrange_debug(ndof, metperm, need_Lg, tree, iinfo )
  If (iinfo < 0 ) Goto 9999
#endif
  
  If( associated(need_Lg)) Deallocate(need_Lg)

  !-----------------------------------------------------------------------------
  ! [4] Compute nnzs for each domain, interface, interior
  !-----------------------------------------------------------------------------
  
  istep = 4
  Call XMPH_Part_compute_domains_nnz(tree, sm_global_A, all_domains, iinfo )
  If (iinfo < 0 ) Goto 9999
  ttmp  = MPI_Wtime() - ttmp
  
  !-----------------------------------------------------------------------------
  ! [5] Save results, free memory, etc.
  !-----------------------------------------------------------------------------
 
  !-------------------------------------------------------------------
  ! [5.1] destroy "tree" ( transmit or free its fields )
  !-------------------------------------------------------------------
  
  ! Transmit & nullify
  all_domains%totinterface = tree%totinterface

  all_domains%metperm     => tree%metperm
  all_domains%domintdof   => tree%domintdof
  all_domains%domstptr    => tree%domstptr

  tree%nbsep        = -1
  tree%totinterface = -1
  tree%toplevel     = -1
  tree%totinterface = -1

  Nullify(tree%metperm)
  Nullify(tree%domintdof)
  Nullify(tree%domstptr)

  ! Free
  Deallocate(tree%node)
  Deallocate(tree%metiperm)
  Deallocate(tree%domptr)

  !-------------------------------------------------------------------
  ! [5.2] Save others
  !-------------------------------------------------------------------

  all_domains%domLg       => domLg
  ana_timing(2) = ttmp
  ana_timing(3) = tsort
  
  !-------------------------------------------------------------------
  ! [5.3] Print error/warnings
  !-------------------------------------------------------------------
9999 Continue
  If ( iinfo /=  0 ) Then
     
     If ( iinfo > 0) msg_class = MSG_WARNING
     If ( iinfo < 0) msg_class = MSG_ERROR
     
     Select Case(istep) 
     Case(1); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
          " while creating the domains")
     Case(21); Call mph_logWithInfo (msg_class,nnza,Trim(rname)//&
          " while allocating memory for the lagrange treatment, size = ")
     Case(22); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
          " while sorting global matrix to CSR")
     Case(31); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
          " while computing the lagrange treatment")
     Case(32); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
          " while exporting the lagrange treatment datas")
     Case(4); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
          " while computing the domains nnz")
     End Select

  End If

  !-------------------------------------------------------------------
  ! [5.4] Set return code
  !-------------------------------------------------------------------

  If ( iinfo == 0 ) info =  0
  If ( iinfo <  0 ) info = -istep
  If ( iinfo >  0 ) info = +istep

End Subroutine XMPH_Part_Build_domains

! [+] routine : XMPH_Part_Split_interface -------------------------------------------
!
!> compute the list of nodes on the interfaces 
!! (their identifier and their associated domain)
!!
!! To construct each domain, we need to identify the interfaces 
!! by analyzing the binary tree.
!! To each node "i" on an interface, 
!! we associate a unique identifier and its domain id,
!! which are respectively stored into 'intrindices(i)' and 'intrfproc(i)' 
!! (range for i is 1..combivtxpcsz).
!!
!! @param[in,out] graph   The matrix graph (on output, unused fields are freed)
!! @param[in,out] tree    The paritioning tree (on output, statistics are appended)
!! @param[   out] domains The domain decomposition description.
!! @param[   out] info    The routine status.
!! 
!! @author Azzam Haidar
!! @author Yohan Lee-tin-yien
!!
!! @todo check the result of allocates
!! @todo implement istep + one goto instead of several gotos.
!! 
Subroutine XMPH_Part_Split_interface & ! intents
     ( graph,                & ! inout
     & tree,                          & 
     & domains,                       & ! out
     & info)

  !* Modules *!


  Use MPH_part_type
  Use mph_log_mod
#if MAPHYS_DEBUG
  Use mph_dbg_mod
#endif
  Implicit None
  
  !* Arguments *!

  Type(maphys_matrix_graph_t)  , intent(inout) :: graph ! deallocate unused fields
  Type(maphys_binary_tree_t) , intent(inout) :: tree ! add statistics
  Type(maphys_domains_t) , intent(out) :: domains
  Integer  , intent(out) :: info  ! subroutine status (<0 = ERROR, 0=SUCCESS, >0 = WARNING)

  !* Local Variables *!

  !-- aliases of structure's fields
  !---- graph fields
  integer                         :: ndof
  integer, dimension (:), pointer :: xadj,adjncy,vwgt
  integer                         :: metalgo

  integer :: maxadj        ! ? a maximum value related to xadj(:) ? - output

  !---- binary tree
  integer                             :: toplevel
  integer                             :: nbsep
  type(maphys_binary_node_t), pointer :: azz (:)
  integer :: totinterface  ! number of vertexes in the separators - output

  integer            , pointer  :: metperm  (:)
  integer            , pointer  :: metiperm (:)

  !---- domains
  integer :: nbdom       ! number of domains     (=nbsep + 1)
  ! integer :: nproc       ! number of processus (=nbdom)

  integer, pointer  :: domintrfdof (:) ! degree of freedom of each domain's interface


  integer :: combivtxpcsz  !< sum on all interface of "myndofintrf" 
  integer, pointer  :: intrfindices(:) ! indexes of interfaces vertex   (one per node per interface)
  integer, pointer  :: intrfproc   (:) ! domain id of interfaces vertex (one per node per interface)

  integer, pointer  :: vtxweight   (:) ! weight (? which one ?) on vertex (size totinterface)
  integer, pointer ::  vtxpos      (:) ! ? vertex position ?              (size totinterface+1)

  integer :: maxprocintrf  ! maximum of interfaces' degrees of freedom 
  integer :: minprocintrf  ! minimum of interfaces' degrees of freedom 

  !-- local variables
  integer :: nproc

  !---- internal error checking
  integer  :: ierr1

  !---- counters
  integer :: i, j, k, l, m, p ! counters      
  integer :: ptt              ! counter on nodes 
  integer :: pt               ! another counter on nodes

  integer :: vtxcnt          ! counter on vertexes
  integer :: vtxcntinit      ! saved counter on vertexes ( at the beginning of a loop )

  integer :: niveaulolsmooth ! smoothing level (see doc_ana.inc -> smoothing)
  integer :: maxlolsmooth    ! ? maximum smoothing level at a certain time ?
  integer :: minlolsmooth    ! ? minimum smoothing level at a certain time ?

  integer :: st              ! starting index
  integer :: ed              ! ending   index

  !---- selected elements
  integer :: procid         ! selected processus
  integer :: sepid          ! selected separator
  integer :: levelid        ! selected level in the binary tree
  integer :: vtx            ! original vertex index                 (1..ndof)
  integer :: vtxid          ! vertex index after maphys permutation (1..ndof) 

  integer :: adjvtx         ! azz: number of my adjacent node in the original indices (=> adjncy(:))
  integer :: belongto       ! selected domain (=> belongingto (:))

  integer :: previndex      ! previous interface node intex (=> intrfindices(:))
  integer :: thisindex      ! selected interface node index (=> intrfindices(:))

  !----- other
  Integer :: iinfo          ! return status of called routine/function
  Integer :: istep          ! flag to a step
  Integer :: msg_class      ! class of a message.
  Integer :: msg_nbelem     ! number of elements to print.
  Integer :: verbosity      ! maphys level of verbosity

  integer :: kflag          ! kflag argument to isort, value forced to 2
  integer :: nothing        ! argument to isort
  integer :: weight         ! ? the weight of a vertex ?
  integer :: dosmooth       ! flag to say if we need to "smooth" (0= no, 1=yes)

  ! colorization the graph's vertexes (= rows (columns) of the global matrix)
  ! if interior, value = processus (=dom) identifier
  ! if on an interface, value = - separator identifier
  integer, allocatable :: belongingto(:)  ! interior or interface (size ndof) 

  integer, allocatable :: vtxrecordisp(:) ! recorded displacement on vertex  (size totinterface)
  integer, allocatable :: existingproc(:) ! related to flowers               (size nproc)
                                          ! (1 = processus exists 0 = do not)   

  ! temporary array to store the processus (= domain) id associated to a vertex
  INTEGER :: vtxpcsz                       !< the size of the temporary array (? estimation ? )
  INTEGER, Pointer :: tempvtxproc (:)  !< the temporary array             (size vtxpcsz)

  Integer, Parameter :: MSGLEN= 2048
  Character(len=MSGLEN) :: msg ! a message.
  

  !- End of header -------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! [1] Initialize local variables
  !-----------------------------------------------------------------------------

  ! 
  Nullify(tempvtxproc)
  
  !-----------------------------------------------------------------------------
  ! [1.1] Init several scalars
  !-----------------------------------------------------------------------------

  ! graph
  ndof    = graph%ndof
  metalgo =  graph%algo
  maxadj  = graph%maxadj
  xadj    => graph%xadj
  adjncy  => graph%adjncy
  vwgt    => graph%vwgt

  ! binary tree
  toplevel = tree%toplevel
  nbsep    = tree%nbsep
  azz      => tree%node
  metperm  => tree%metperm
  metiperm => tree%metiperm

  ! local variable
  nbdom = nbsep + 1
  nproc = nbdom

  info = 0            
  niveaulolsmooth = 0 

  !-----------------------------------------------------------------------------
  ! [1.2] Compute belongingto (links each row/column to a domain or a separator)
  !-----------------------------------------------------------------------------
  !
  ! For i=1,ndof
  ! If belongingto(i) > 0 means row/col "i" belongs to domain    " belongingto(i)" 
  ! If belongingto(i) < 0 means row/col "i" belongs to separator "-belongingto(i)" 
  !
  
  ! Allocate
  istep = 121
  ALLOCATE(belongingto(ndof),STAT=ierr1)
  Do i=1, ndof
     belongingto(i)=999
  End Do

  ! Compute
  istep = 122
  totinterface = 0
  do i=1, nbsep 
     totinterface=totinterface+(azz(i)%seped-azz(i)%sepst+1)
     if(azz(i)%rson.eq.-1) then 
        procid = azz(i)%rdomid 
        do j=azz(i)%rgst,azz(i)%rged
           belongingto(metperm(j)) = procid
        enddo
     endif
     if(azz(i)%lson.eq.-1) then 
        procid = azz(i)%ldomid 
        do j=azz(i)%lgst,azz(i)%lged
           belongingto(metperm(j)) = procid
        enddo
     endif
     sepid = -i 
     do j=azz(i)%sepst,azz(i)%seped
        belongingto(metperm(j)) = sepid
     enddo
  enddo

  ! Print statistics

  Call mph_logWithInfo&
       (MSG_VERBOSE,totinterface, "Domain Decomposition -> total size of the interface =")
  Call mph_logWithInfo&
       (MSG_VERBOSE,maxadj,  "Domain Decomposition -> maximum of adjacency =")

#if MAPHYS_DEBUG  
  ! Check
  istep = 123
  do i=1, ndof 
     if(belongingto(i).ge.nproc)Then
        iinfo = -1
        Goto 9999
     End if
  enddo
#endif

#if MAPHYS_DEBUG
  ! Write separator to a file
  istep = 124
  Call MPH_dbg_init
  Call MPH_dbg_set_file("separator.mtx")

  write(UNIT=dbg_unit,FMT='(A)')"%%MatrixMarket matrix coordinate REAL general" 
  write(UNIT=dbg_unit,FMT='(I16,I16,I16)') totinterface,1,totinterface
  k=0
  do i=1, nbsep 
     do j=azz(i)%sepst,azz(i)%seped
        k=k+1
        write(unit=dbg_unit,FMT='(I16,I16,I16)') k,1,metperm(j)
     enddo
  enddo

#endif

  !-----------------------------------------------------------------------------
  ! [1.3] Estimate necessary memory size
  !-----------------------------------------------------------------------------
  !> @warning estimations of memory (Yohan)
  !! As Azzam told me those formula are purely empirical  
  !!
  combivtxpcsz =  4 * totinterface  ! this array maybe reallocated
  vtxpcsz  = 36 * (maxadj+1)        ! this array maybe reallocated. 

  !-----------------------------------------------------------------------------
  ! [1.4] Allocate arrays
  !-----------------------------------------------------------------------------

  istep = 14

  Allocate(intrfindices(combivtxpcsz),STAT=iinfo)
  If (iinfo > 0) iinfo = -iinfo
  If (iinfo < 0) Goto 9999

  Allocate(intrfproc(combivtxpcsz)   ,STAT=iinfo)
  If (iinfo > 0) iinfo = -iinfo
  If (iinfo < 0) Goto 9999

  Allocate(tempvtxproc(vtxpcsz)  ,STAT=iinfo)
  If (iinfo > 0) iinfo = -iinfo
  If (iinfo < 0) Goto 9999

  Allocate(vtxpos(totinterface+1),STAT=iinfo)
  If (iinfo > 0) iinfo = -iinfo
  If (iinfo < 0) Goto 9999

  Allocate(vtxweight(totinterface),STAT=iinfo)
  If (iinfo > 0) iinfo = -iinfo
  If (iinfo < 0) Goto 9999

  Allocate(vtxrecordisp(totinterface),STAT=iinfo)
  If (iinfo > 0) iinfo = -iinfo
  If (iinfo < 0) Goto 9999

  Allocate(domintrfdof(nproc),STAT=iinfo)
  If (iinfo > 0) iinfo = -iinfo
  If (iinfo < 0) Goto 9999

  Allocate(existingproc(nproc),STAT=iinfo)
  If (iinfo > 0) iinfo = -iinfo
  If (iinfo < 0) Goto 9999


  !-----------------------------------------------------------------------------
  ! [2] Form two vectors [intrfindices intrfproc]
  !-----------------------------------------------------------------------------

  ! ------------------------------------
  !       SEPARATOR ASSOCIATION
  ! ------------------------------------
  !      GO THROUGH SEPARATOR TREE 
  !  from level 1 to toplevel and assign 
  !  each vertex to the subdomain that share.
  !  So construct two vectors [intrfindices intrfproc]
  !  where: each ith value correspond to a 
  !  couple of value that are (an interface 
  !  index and its corresponding processor or subdomain)
  !  note that here the interface indices are 
  !  repeted according to there belonging processors
  ! ------------------------------------

  !-----------------------------------------------------------------------------
  ! [2.1] Init arrays
  !-----------------------------------------------------------------------------

  Do i=1, combivtxpcsz
     intrfindices(i)=0
  End Do
  
  Do i=1, combivtxpcsz    
     intrfproc(i)=0
  End Do

  ptt       =  0
  vtxcnt    =  0
  Do i=1, totinterface + 1
     vtxpos(i) = -1
  End Do
  vtxpos(1) =  1

  Do i=1, nproc
     domintrfdof(i)=  0
  End Do

  Do i=1, totinterface
     vtxweight(i)=0
  End Do

  
  ! I will use it to find the displacment of each vtx in vtxpos.                                                
  ! vtxpos is like iacsr for the starting point of each vtx in [intrfindices intrfproc].                        
  
  ! BECAUSE the treatment of vertex is not in order (1 2 3) so vtxpos(i)                                        
  ! give the starting pointer to the ith treated vtx and not the ith interface of [Abb].                        
  ! So I use vtxrecordisp to redirect the vtx index to their treatment order to find their corresponding vtxpos.

  Do i=1, totinterface
     vtxrecordisp(i)=0 
  End Do

  !-----------------------------------------------------------------------------
  ! [2.2] Compute the arrays
  !-----------------------------------------------------------------------------

  do levelid=1,toplevel 
     ! ------------------------------------
     !        FOR EACH SEPARATOR
     ! of this levelid find corresponding domainid
     ! ------------------------------------
     do j=1,nbsep
        if(azz(j)%level.ne.levelid) goto 320  ! not the level in question so omit it  !lty: == CYCLE
        ! ------------------------------------
        ! Azzam: Modification 05/01/2010
        !  FIRST LEVEL SEPARATOR: SO ASSOCIATE 
        !   IT TO THE CORRESPONDING SUBDOMAIN 
        !         THAT IT SPLIT THEM
        ! ------------------------------------
        if(azz(j)%level.eq.1)then 
           ! -------------------------------------------------------------
           ! to keep order in proc assignment from smallest to biggest, 
           ! choose between storing the right or the left subdomain first
           ! -------------------------------------------------------------
           if(azz(j)%rdomid.lt.azz(j)%ldomid)then 
              do k=azz(j)%sepst,azz(j)%seped
                 ! associate right subdomain
                 ptt               = ptt+1
                 if(ptt.gt.combivtxpcsz)Then
                    Call combivtxpcsz_realloc (&
                         combivtxpcsz, ptt, intrfindices,intrfproc, info)
                    MPH_ONFAILURE_GOTO9999(info)
                 End if
                 procid            = azz(j)%rdomid 
                 intrfindices(ptt) = k 
                 intrfproc(ptt)    = procid
                 domintrfdof(procid+1) = domintrfdof(procid+1)+1
                 ! associate left subdomain
                 ptt               = ptt+1
                 if(ptt.gt.combivtxpcsz)Then
                    Call combivtxpcsz_realloc (&
                         combivtxpcsz, ptt, intrfindices,intrfproc, info)
                    MPH_ONFAILURE_GOTO9999(info)
                 End if
                 procid            = azz(j)%ldomid 
                 intrfindices(ptt) = k 
                 intrfproc(ptt)    = procid
                 domintrfdof(procid+1) = domintrfdof(procid+1)+1
                 ! I have associate 2 subdomain to this vtx
                 vtxcnt            = vtxcnt +1
                 vtxrecordisp(k)   = vtxcnt
                 vtxpos(vtxcnt+1)  = ptt+1 
                 vtxweight(k)      = 2
              enddo ! END DO OF FOR EACH VERTEX (==>ALL VERTEX OF THIS SEPARATOR ARE DONE)
           else 
              do k=azz(j)%sepst,azz(j)%seped
                 ! associate left subdomain
                 ptt               = ptt+1
                 if(ptt.gt.combivtxpcsz)Then
                    Call combivtxpcsz_realloc (&
                         combivtxpcsz, ptt, intrfindices,intrfproc, info)
                    MPH_ONFAILURE_GOTO9999(info)
                 End if
                 procid            = azz(j)%ldomid 
                 intrfindices(ptt) = k 
                 intrfproc(ptt)    = procid
                 domintrfdof(procid+1) = domintrfdof(procid+1)+1
                 ! associate right subdomain
                 ptt               = ptt+1
                 if(ptt.gt.combivtxpcsz)Then
                    Call combivtxpcsz_realloc (&
                         combivtxpcsz, ptt, intrfindices,intrfproc, info)
                    MPH_ONFAILURE_GOTO9999(info)
                 End if
                 procid            = azz(j)%rdomid 
                 intrfindices(ptt) = k 
                 intrfproc(ptt)    = procid
                 domintrfdof(procid+1) = domintrfdof(procid+1)+1
                 ! I have associate 2 subdomain to this vtx
                 vtxcnt            = vtxcnt +1
                 vtxrecordisp(k)   = vtxcnt
                 vtxpos(vtxcnt+1)  = ptt+1 
                 vtxweight(k)      = 2
              enddo ! END DO OF FOR EACH VERTEX (==>ALL VERTEX OF THIS SEPARATOR ARE DONE)
           endif ! endif of right or left subdomain
           ! -------------------------------------------------------------
           ! -------------------------------------------------------------
           goto 320 !lty: cycle
        endif  ! ENDIF FOR FIRST LEVEL SEPARATOR
        ! ------------------------------------
        ! END Modification 05/01/2010
        ! ------------------------------------

        !            print*,'treatment of separator number',j
        ! ------------------------------------
        !  FOR EACH VERTEX OF THIS SEPARATOR
        ! ------------------------------------
        ! for each vertex (vtx) of this separator  
        ! find its adjacency and associate this  
        ! vtx to the corresponding processors
        ! ------------------------------------
        st=ptt+1
        vtxcntinit=vtxcnt
        minlolsmooth=0
        do k=azz(j)%sepst,azz(j)%seped

           pt=0
           vtx = metperm(k) ! number of this vertex in the original indices 
           ! ------------------------------------
           !   FOR EACH ADJCENT OF THIS VERTEX
           ! ------------------------------------
           ! FOR EACH ADJCENT OF THIS VERTEX:
           ! IF belongto>0 ==> adjvtx is an interior node ==> so I am interface for the corresponding subdomain
           ! ELSE adjvtx belong to a separator  
           !      IF the level of this separator is less than my level associate to myself all subdomains that this adjvtx belong to.
           !      ELSE IF the level is bigger than my level, keep it to my father separator at next round.
           !      ELSE IF the level is equal to my level, I should take all subdomains that this adjvtx belong to. 
           !              BUT I don't know if this adjvtx is already done. 
           !              So if this adjvtx is done ==> its easy take all subdomains that it belong to.
           !              else I have to wait that I finish this separator then make a round on this separator to take this adjvtx into account.
           ! --------------
           !   NIVEAU -1-
           ! --------------
           do l=xadj(vtx),xadj(vtx+1)-1  
              adjvtx   = adjncy(l)           ! number of my adjacent node in the original indices 
              belongto = belongingto(adjvtx) ! number of domain or separator that my adjacent node (adjvtx) belong to
              if(belongto.ge.0) then         ! I am interface for the corresponding subdomain
                 pt=pt+1                     
                 If(pt.gt.vtxpcsz)Then
                    Call vtxpcsz_realloc (&
                         vtxpcsz, pt, tempvtxproc, info)
                    MPH_ONFAILURE_GOTO9999(info)
                 End If
                 tempvtxproc(pt) = belongto   
              elseif(azz(-belongto)%level.lt.levelid) then ! adjvtx belong to a separator less than my level 
                 vtxid=metiperm(adjvtx)                    ! index of adjvtx in the new perm [separator interior]
                 vtxid=vtxrecordisp(vtxid)                 ! redirection to indiquate the position of vtx in the order of treatment.
                 do m=vtxpos(vtxid),vtxpos(vtxid+1)-1
                    pt=pt+1
                    If(pt.gt.vtxpcsz)Then
                       Call vtxpcsz_realloc (&
                            vtxpcsz, pt, tempvtxproc, info)
                       MPH_ONFAILURE_GOTO9999(info)
                    End If
                    tempvtxproc(pt) = intrfproc(m) 
                 enddo
              elseif(azz(-belongto)%level.ge.levelid) then 
                 vtxweight(k)     = -1 ! put -1 to mean that this vtx need niveau 2 
              endif
           enddo ! END DO OF FOR EACH ADJCENT (==>ALL ADJCENT VERTEX ARE DONE)
!lty: stop

           ! Now tempvtxproc contain all the processor that are associated with vtx (with duplication)
           ! so associate this vtx to its corresponding processors by order without duplication 
           !               if(pt.eq.0) goto 978 ! impossible. vertex is isolated and does not belong to any domain
           If(pt.gt.vtxpcsz)Then
              Call vtxpcsz_realloc (&
                   vtxpcsz, pt, tempvtxproc, info)
              MPH_ONFAILURE_GOTO9999(info)
           End If
           !if(ptt+pt.gt.combivtxpcsz) goto 976 ! 
           kflag = 1
           CALL ISORT(tempvtxproc,nothing,pt,kflag)
           ptt               = ptt+1
           weight            = ptt
           If(ptt.gt.combivtxpcsz)Then
              Call combivtxpcsz_realloc (&
                   combivtxpcsz,ptt, intrfindices,intrfproc, info)
              MPH_ONFAILURE_GOTO9999(info)
           End if
           procid            = tempvtxproc(1)
           intrfindices(ptt) = k 
           intrfproc(ptt)    = procid
           domintrfdof(procid+1) = domintrfdof(procid+1) + 1
           do p=2,pt
              if(tempvtxproc(p).ne.tempvtxproc(p-1))then
                 ptt               = ptt+1
                 If(ptt.gt.combivtxpcsz)Then
                    Call combivtxpcsz_realloc (&
                         combivtxpcsz, ptt, intrfindices,intrfproc, info)
                    MPH_ONFAILURE_GOTO9999(info)
                 End if
                 procid            = tempvtxproc(p)
                 intrfindices(ptt) = k 
                 intrfproc(ptt)    = procid
                 domintrfdof(procid+1) = domintrfdof(procid+1)+1
              endif
           enddo
           vtxcnt            = vtxcnt +1
           vtxrecordisp(k)   = vtxcnt
           vtxpos(vtxcnt+1)  = ptt+1 
           weight            = ptt-weight+1
           if(vtxweight(k).eq.0)then
              vtxweight(k)= weight
           else
              vtxweight(k)= -weight
           endif
           !              ----------------------------------
           !              modif pour lancer niveau 2 lorsque le weight est < que 2
           if(weight.lt.2)minlolsmooth=1
           !              ----------------------------------
        enddo ! END DO OF FOR EACH VERTEX (==>ALL VERTEX OF THIS SEPARATOR ARE DONE)

        l=azz(j)%seped-azz(j)%sepst+1
        m=vtxcnt-vtxcntinit
        if(m.ne.l) goto 957

        !           ----------------------------------
        !           ----------------------------------
        !            NIVEAU -2-   LISSAGE-SMOOTHING
        !           ----------------------------------
        !           ----------------------------------
        !           Once finishing the first round of this separator vtx's,
        !           make another round (smoothing) to find the extra couple (vtx proc).
        !           Once finishing this, sort intrfindices to keep order in vtx.
        !           then update the database coz it is used by the next level of separator. 
        !           ----------------------------------
        !           Find extra proc
        niveaulolsmooth=0
        maxlolsmooth=10
        minlolsmooth=1
        dosmooth=minlolsmooth
        DO WHILE &
             ((dosmooth.NE.0).AND.(niveaulolsmooth.LT.maxlolsmooth)&
             .AND.(minlolsmooth.GE.1))
           minlolsmooth=minlolsmooth-1
           if(minlolsmooth.le.0)dosmooth=0
           ed=ptt

           do k=azz(j)%sepst,azz(j)%seped
              if(vtxweight(k).lt.0)then ! vtx couple association is not finished at niveau -1- and need niveau-2- for smoothing 
                 vtxweight(k)     =  0                 ! put 0 to mean that niveau 2 on this vtx is done
                 !                    find the existing proc that this vtx belong to them 
                 Do l=1, nproc
                    existingproc(l)=0
                 End Do

                 vtxid=vtxrecordisp(k)                 
                 do l=vtxpos(vtxid),vtxpos(vtxid+1)-1
                    existingproc(intrfproc(l)+1)=1
                 enddo
                 !                    start smoothing
                 vtx = metperm(k) ! number of this vertex in the original indices 
                 do l=xadj(vtx),xadj(vtx+1)-1  
                    adjvtx   = adjncy(l)           ! number of my adjacent node in the original indices 
                    belongto = belongingto(adjvtx) ! number of domain or separator that my adjacent node (adjvtx) belong to
                    if(belongto.lt.0) then         ! this is an interface adjncy that belong to a separator.
                       if(azz(-belongto)%level.eq.levelid) then ! adjncy belong to my separator ==> do the niveau 2 of smoothing
                          vtxid=metiperm(adjvtx)                ! index of adjvtx in the new perm [separator interior]
                          vtxid=vtxrecordisp(vtxid)             
                          do m=vtxpos(vtxid),vtxpos(vtxid+1)-1
                             procid=intrfproc(m)
                             if(existingproc(procid+1).eq.0) then   
                                ptt               = ptt+1
                                If(ptt.gt.combivtxpcsz)Then
                                   Call combivtxpcsz_realloc (&
                                        combivtxpcsz, ptt, intrfindices,intrfproc, info)
                                   MPH_ONFAILURE_GOTO9999(info)
                                End if
                                intrfindices(ptt) = k 
                                intrfproc(ptt)    = procid
                                existingproc(procid+1)=1
                                domintrfdof(procid+1) = &
                                     domintrfdof(procid+1)+1
                             endif
                          enddo
                       elseif(azz(-belongto)%level.gt.levelid) then 
                          vtxweight(k)     = -1 
                       endif
                    endif
                 enddo
              endif
           enddo ! END DO OF FOR EACH VERTEX (==>ALL VERTEX OF THIS SEPARATOR ARE DONE niveau 2) 


           ! IF    I add some couple (vtx proc) to the list of this separator, I need to sort and update
           ! ELSE  I need just to update weight from negatif to positif value
!lty: stop

           if(ptt.gt.ed)then ! it mean that I add some couple (vtx proc) to the list (intrfindices,intrfproc) of this separator.
              ! order them by vertex and update the database
              kflag = 2
              p=ptt-st+1     ! length of couple (vtx proc) of this separator to be ordered.
              ! p should be greather than 1 coz a separator is formed by at least one vtx and a vtx should belong to at least 2 proc
              CALL ISORT(intrfindices(st),intrfproc(st),p,kflag)

              l=st
              pt=st
              previndex = intrfindices(st)  ! index of the vtx in the new perm [separator interior]
              weight    = 1
              vtxcnt    = vtxcntinit
              do l=st+1,ptt
                 thisindex = intrfindices(l)
                 if(thisindex.eq.previndex)then             ! old vtx increment counter
                    pt                = pt+1
                    weight            = weight+1
                 else                                       ! new vtx update, then initialize counter and previndex
                    if(vtxweight(previndex).eq.0)then
                       vtxweight(previndex)      = weight
                    elseif(vtxweight(previndex).eq.-1)then
                       vtxweight(previndex)      = -weight
                    endif
                    if(weight.eq.1) then                    ! impossible. vertex  belong to 1 domain so how can be an interface?
                       dosmooth=1
                       vtxweight(previndex)=-1
                    endif
                    vtxcnt                    = vtxcnt +1
                    vtxrecordisp(previndex)   = vtxcnt
                    vtxpos(vtxcnt+1)          = pt+1 
                    ! initialize counter and previndex
                    previndex                 = thisindex
                    pt                        = pt+1
                    weight                    = 1
                 endif
              enddo
              ! take last couple (vtx proc) into consideration coz the last vtx should 
              ! be the same as the previous one so it is not updated for this vtx
              thisindex = intrfindices(ptt)
              if(l-1.eq.ptt) then                           ! last value in the list
                 if(thisindex.eq.previndex)then             ! the last vtx is equal to the one before so just update
                    if(vtxweight(previndex).eq.0)then
                       vtxweight(previndex)      = weight
                    elseif(vtxweight(previndex).eq.-1)then
                       vtxweight(previndex)      = -weight
                    endif
                    if(weight.eq.1) then                    ! impossible. vertex  belong to 1 domain so how can be an interface?
                       dosmooth=1
                       vtxweight(previndex)=-1
                    endif
                    vtxcnt                    = vtxcnt +1
                    vtxrecordisp(previndex)   = vtxcnt
                    vtxpos(vtxcnt+1)          = pt+1 
                 else ! mean that this last vtx is alone and has been done in the do loop 
                    dosmooth=1
                    vtxweight(previndex)=-1
                 endif
              else ! IMPOSSIBLE COZ at output from the do loop l=ptt+1
                 Call mph_log(MSG_ERROR,'IMPOSSIBLE FALSE DO LOOP')
                 goto 555
              endif
              if(pt.ne.ptt)goto 959
           else ! it mean that I doesn't add any couple, so just update the weight of this separator

              do k=azz(j)%sepst,azz(j)%seped
                 vtxid=vtxrecordisp(k)                 ! redirection to indiquate the position of vtx in the order of treatment.
                 vtxweight(k)=0
                 if(vtxweight(k).eq.0)then
                    do l=vtxpos(vtxid),vtxpos(vtxid+1)-1
                       vtxweight(k)=vtxweight(k)+1
                    enddo
                 elseif(vtxweight(k).eq.-1)then
                    do l=vtxpos(vtxid),vtxpos(vtxid+1)-1
                       vtxweight(k)=vtxweight(k)-1
                    enddo
                 endif
                 if(vtxweight(k).eq.1)goto 951
              enddo
              dosmooth=0

              Write(msg,'(A,I5,A,I10,A)') 'for separator Id',j,&
                   ': smooth level 2 stopped at level',&
                   niveaulolsmooth+1,' due to a empty smoothed loop '
              Call mph_log(MSG_WARNING,msg)

           endif
           l=azz(j)%seped-azz(j)%sepst+1
           m=vtxcnt-vtxcntinit
           if(m.ne.l) goto 957

           niveaulolsmooth=niveaulolsmooth+1
           if((niveaulolsmooth.ge.maxlolsmooth).and.(dosmooth.ne.0))Then
              Write(msg,'(A,I5,A)') &
                   'for separator Id',j, 'smooth level 2 stopped at maxsmooth'
              Call mph_log(MSG_WARNING,msg)
           End if

        ENDDO
        !           ---------------------------------
        !           END NIVEAU -2-
        !           ---------------------------------

320     continue 

     enddo ! END DO OF FOR EACH SEPARATOR (==>ALL SEPARATOR of THE LEVEL IN QUESTION ARE DONE)
  enddo ! END FOR ALL LEVEL

  Call mph_logWithInfo(MSG_VERBOSE,&
       niveaulolsmooth,'Domain Decomposition -> 2nd smoothing -> final level =')

  do k=1,totinterface
     if(vtxweight(k).lt.0)then
        vtxweight(k)      = -vtxweight(k)
     endif
     if(vtxweight(k).eq.1)goto 950
     if(vtxweight(k).eq.0)goto 949
  enddo
  combivtxpcsz=ptt
  minprocintrf=minval(domintrfdof)
  maxprocintrf=maxval(domintrfdof)
  if(vtxcnt.ne.totinterface) goto 958

  !*     !--------------------
  !*     ! form back metiperm  NO NEED TO THIS COZ IT IS NOT ANYMORE MODIFIED ABOVE IN THIS VERSION
  !*     !--------------------
  !      do k=1,totinterface
  !         vtx = metperm(k) ! number of this vertex in the original indices 
  !         metiperm(vtx)=k
  !      enddo
  !*     !--------------------


  !-----------------------------------------------------------------------------
  ! [ ] Print the statistics on local interface
  !-----------------------------------------------------------------------------

  ! Give summary statistics

  Call mph_logWithInfo(MSG_STD,maxval(domintrfdof),&
  'Domain Decomposition -> local interface -> maximum size = ')

  Call mph_logWithInfo(MSG_STD,minval(domintrfdof),&
  'Domain Decomposition -> local interface -> minimum size =')

  Call mph_logWithInfo(MSG_STD,minval(vtxweight),&
  'Domain Decomposition -> local interface -> maximum weight =')

  ! Give several values of domintrdof
  Call mph_log_GetVerbosity(verbosity)

  msg_nbelem = 0
  If ( verbosity == MSG_STD     ) msg_nbelem = min(10,nproc)
  If ( verbosity >= MSG_VERBOSE ) msg_nbelem = nproc
  If ( verbosity == MSG_STD     ) msg_class  = MSG_STD
  If ( verbosity >= MSG_VERBOSE ) msg_class  = MSG_VERBOSE

  If ( msg_nbelem > 0 ) Then
     Do i=1,msg_nbelem
        Write( msg,'(2A,I5,A,I10)') 'Domain Decomposition ->',&
             ' local interface [',i-1,'] -> size = ',&
             domintrfdof(i)
        Call mph_log(msg_class, msg)
     End Do
  End If

  If ( msg_nbelem > 0 ) Then
     Do i=1,msg_nbelem
        Write( msg,'(2A,I5,A,I10)') 'Domain Decomposition ->',&
             ' local interior [',i-1,'] -> size = ',&
             tree%domintdof(i)
        Call mph_log(msg_class, msg)
     End Do
  End If

  
  !
  !
  !     =========================================
  !      finish forming [intrfindices intrfproc]
  !     =========================================
  !      ATTENTION [intrfindices intrfproc] are sorted according
  !      to the indices. intrfindices is sorted and carry intrfproc along 
  !     =========================================
  
  
  IF(metalgo.eq.3) DEALLOCATE(vwgt)
  DEALLOCATE(xadj) 
!      DEALLOCATE(adjncy)
!      DEALLOCATE(belongingto)

  DEALLOCATE(tempvtxproc)
  DEALLOCATE(vtxrecordisp)
  DEALLOCATE(existingproc)


  ! associate local variables to output structures
  ! ==============================================
  tree%totinterface =  totinterface 

  domains%nbdom        = nbdom       
  domains%maxprocintrf = maxprocintrf
  domains%minprocintrf = minprocintrf

  domains%combivtxpcsz =  combivtxpcsz 

  domains%intrfindices=>  intrfindices 
  domains%intrfproc   =>  intrfproc    

  domains%domintrfdof =>  domintrfdof  
  domains%vtxweight   =>  vtxweight    
  domains%vtxpos      =>  vtxpos       

#if MAPHYS_DEBUG
  Call MPH_dbg_init

  Call MPH_dbg_set_file("domains.txt")
  Write(dbg_unit,*) "totinterface =", domains%totinterface 
  Write(dbg_unit,*) "nbdom        =", domains%nbdom
  Write(dbg_unit,*) "maxprocintrf =", domains%maxprocintrf
  Write(dbg_unit,*) "minprocintrf =", domains%minprocintrf
  Write(dbg_unit,*) "combivtxpcsz =", domains%combivtxpcsz
  
  Write(dbg_unit,*) "intrfindices =", domains%intrfindices(1:combivtxpcsz)
  Write(dbg_unit,*) "intrfproc    =", domains%intrfproc   (1:combivtxpcsz)

  Write(dbg_unit,*) "domintrfdof  =", domains%domintrfdof (1:nproc)
  Write(dbg_unit,*) "vtxweight    =", domains%vtxweight   (1:totinterface)
  Write(dbg_unit,*) "vtxpos       =", domains%vtxpos      (1:totinterface+1)
#endif

  !---------------------------------------------------------------------------
  ! [*.*] Error handling
  !---------------------------------------------------------------------------

9999 Continue

    ! Print error messages

    If (iinfo /= 0) Then
       If ( iinfo > 0) msg_class = MSG_WARNING
       If ( iinfo < 0) msg_class = MSG_ERROR
     
       Select Case(istep) 
       Case(123); Write(msg,'(A,I20,I20)')                                   &
            'bad belongingto (index value) = ',i,belongingto(i)
       Case(14 ); Write(msg,'(A)')                                           &
            'failed to allocate components of "domains"'

#if MAPHYS_DEBUG
       Case(978); Write(msg,'(A,A,I10)')                                     &  
            ' ERROR 978: the number of subdomain that a separator vertex',   & 
            ' belong to is equal to zero ',vtx
#endif
       Case(959); write(msg,FMT='(3A,2I10)')                                 &
            ' Azz_const_struct ERROR 959 after updating at niveau -2- the ', &
            'value of the added couple pt after sorting is different from ', &
            'the value of added ptt at smoothing.',pt,ptt
       Case(958); write(msg,FMT='(A60,I10,I10)')                             &
            'Azz_const_struct ERROR 958 (vtxcnt.ne.totinterface)',           &
            vtxcnt,totinterface
       Case(957); write(msg,FMT='(A,I10,A,2I10,1/)')                         &
            'Azz_const_struct ERROR 957 vtxcnt of separator ',j,             &
            'not equal to the number of node in the separator',m,l
       Case(951); write(msg,'(A,A,I16)')                                     &
            ' ERROR 951: the number of subdomain that a separator vertex',   &
            ' belong to is equal to one in the smoothed niveau2',k
       Case(950); write(msg,'(A,A,I16)')                                     &
            ' ERROR 950: the number of subdomain that a separator vertex',   &
            ' belong to is equal to one',k
       Case(949); write(msg,'(A,I16)')                                       &
            ' ERROR 949: vtxweight(k)= 0 at niveau2bis',k 
       Case(555); Write(msg,'(A)'    )                                       &
            'Internal error = 555 '
       Case Default; Write(msg,'(A,I16)')                                    &
            'Unreferenced internal error =', istep
       End Select

       Call mph_log(msg_class, Trim(msg))

    End If

    ! Set return code

    If ( iinfo == 0 ) info = 0
    If ( iinfo <  0 ) info = -istep
    If ( iinfo >  0 ) info = +istep

    ! Returning
    
    Return

    ! binding to errors
#if MAPHYS_DEBUG
978 istep=978 ; iinfo = -1; Goto 9999
#endif

959 istep=959 ; iinfo = -1; Goto 9999
958 istep=958 ; iinfo = -1; Goto 9999
957 istep=957 ; iinfo = -1; Goto 9999
951 istep=951 ; iinfo = -1; Goto 9999
950 istep=950 ; iinfo = -1; Goto 9999
949 istep=949 ; iinfo = -1; Goto 9999
555 istep=555 ; iinfo = -1; Goto 9999

  End Subroutine XMPH_Part_Split_interface

  !> reallocate the array which size depend on combivtxpcsz
  Subroutine combivtxpcsz_realloc &
       ( size, minsize, intrfindices,intrfproc, info)

    !* Module *!
    Use XMPH_sparse_matrix_mod, Only : XMPH_irealloc
    Use MPH_Log_mod
    Implicit None
    
    !* Arguments *!
    MPH_INT, Intent(inout) :: size
    MPH_INT, Intent(in   ) :: minsize
    MPH_INT, Pointer, Intent(inout) :: intrfindices(:), intrfproc(:)
    Integer, Intent(  out) :: info

    !* Local variables *!
    MPH_INT ::  newsize

    ! End of header -----------------------------------------------------------------------

    ! Inform user
    newsize = Max(minsize,Int(1.4 * size + 1,MPH_INTKIND))
    Call mph_logWithInfo(MSG_WARNING,newsize,&
         "ptt bigger than combivtxpcsz. Several arrays are to be reallocated.new combivtxpcsz:") 

    Call XMPH_irealloc(intrfindices,newsize,size,info)
    If (info < 0) Goto 9999
    Call XMPH_irealloc(intrfproc,newsize,size,info)
    If (info < 0) Goto 9999
    size = newsize

    ! Inform user
9999 Continue
    If (info < 0) Then
       Call mph_log(MSG_ERROR,"Not enough memory or invalid size.") 
    Else
       Call mph_log(MSG_WARNING,"Memory reallocated.") 
    End If

  End Subroutine combivtxpcsz_realloc

  !> same as combivtxpcz_realloc but for vtxpcsz
  Subroutine vtxpcsz_realloc &
       (size, minsize, tempvtxproc, info)

    !* Module *!
    Use XMPH_sparse_matrix_mod, Only : XMPH_irealloc
    Use MPH_Log_mod
    Implicit None

    !* Arguments *!
    MPH_INT, Intent(inout) :: size
    MPH_INT, Intent(in   ) :: minsize
    MPH_INT, Pointer, Intent(inout) :: tempvtxproc(:)
    Integer, Intent(  out) :: info

    !* Local variables *!
    MPH_INT ::  newsize

    ! End of header -----------------------------------------------------------------------

    ! Inform user
    newsize = Max(minsize,Int(1.4 * size + 1,MPH_INTKIND))
    Call mph_logWithInfo(MSG_WARNING,newsize,&
         "pt bigger than vtxpcsz. Several arrays are to be reallocated.new vtxpcsz:") 
    Call XMPH_irealloc(tempvtxproc,newsize,size,info)
    If (info < 0) Goto 9999
    size = newsize

    ! Inform user
9999 Continue
    If (info < 0) Then
       Call mph_log(MSG_ERROR,"Not enough memory or invalid size.") 
    Else
       Call mph_log(MSG_WARNING,"Memory reallocated.") 
    End If

  End Subroutine vtxpcsz_realloc

  !--------------------------------------------------------------------------------------

! [+] routine : XMPH_Part_sort_GlobalMatrix_CSR -----------------------------------------
!
!> convert the input global matrix into CSR format and detect which row/columns are Lagrange
!! 
!!
!!
!! @param [in,out] sm_global 
!!        The matrix to be sorted into CSR
!!
!! @param [in,out] need_Lg   
!!        Integer Array (size = order of matrix) allocated by the caller.
!!        - On input, its values are irelevant.
!!        - On output, specifies if the row/column needs Lagrange treatment.
!!         Indeed, if entry i, is on the diagonal and is equal to zero, need_Lg(i) = 0; 
!!         else need_Lg(i) = 1.
!!        .
!! @param [   out] tsort
!!        The time (in seconds) to sort sm_global into CSR.      
!!
!! @param [   out] info 
!!        The routine status
!!
!! @author Azzam Haidar
!! @author Yohan Lee-tin-yien
!!
Subroutine XMPH_Part_sort_GlobalMatrix_CSR( sm_global, need_Lg, tsort, info )

  !* Modules *!
  Use XMPH_sparse_matrix_mod
  Use MPH_part_type
  Use mph_log_mod
  implicit none

  !* Arguments *!
  type(XMPH_sparse_matrix_t) , intent(inout) :: sm_global
  MPH_INT  , pointer , intent(inout) :: need_Lg(:)

  real(kind=8)          , intent(  out) :: tsort
  integer               , intent(  out) :: info
  
  !* External functions *!
  real(kind=8), external :: MPI_Wtime

  !* Local variables *!

  ! Scalars 
  Integer               :: iinfo     ! status of called routines
  Integer               :: istep     ! internal step
  Integer               :: msg_class ! class of messages to print
  Integer               :: kflag     ! Flag for IPSORT
  MPH_INT            :: ndof      ! sm_global%m
  MPH_INT            :: nnza      ! sm_global%nnz
  MPH_INT            :: i, ivp    ! counters

  ! Arrays
  MPH_INT  , pointer :: ia(:)    ! sm_global row indices.
  MPH_INT  , pointer :: ja(:)    ! sm_global column indices
  MPH_INT  , pointer :: iacsr(:) ! sm_global  compressed rows indices
  XMPH_FLOAT, pointer ::  a(:)    ! sm_global values

  ! String
  Character (len=MAPHYS_STRL), Parameter :: rname = "XMPH_Part_SORT_GLOBALMATRIX_CSR"

  ! End of header --------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! [1] Init
  !-----------------------------------------------------------------------------

  ndof =  sm_global%m
  nnza =  sm_global%nnz
  ia   =>  sm_global%i
  ja   =>  sm_global%j
  a    =>  sm_global%v

  Do i=1,nnza
     need_Lg(i)=0
  End Do

  istep = 1
  Allocate(iacsr(ndof+1), STAT=iinfo )
  If (iinfo >  0 ) iinfo = -iinfo
  If (iinfo <  0 ) Goto 9999
  
  !-----------------------------------------------------------------------------
  ! [2] Sort sm_global according to the row indices (in increasing order)
  !-----------------------------------------------------------------------------

  !     ======================================================
  !                 SORT IA JA A and form IACSR
  !
  !     ATTENTION HERE AT THIS STEP (BEFORE SORTING)
  !     IA JA SHOULD BE PERMUTED ACCORDING TO MY PERMUTATION 
  !     that is:        
  !     separator_row [ Abb  Abi1   Abi2    Abi3   Abi4   Abi5]
  !     interior__row [ Aib1 Aii1    0       0      0       0 ]
  !     interior__row [ Aib2   0    Aii2     0      0       0 ]
  !     interior__row [ Aib3   0     0      Aii3    0       0 ]
  !     interior__row [ Aib4   0     0       0     Aii4     0 ]
  !     interior__row [ Aib5   0     0       0      0     Aii5]
  !     ======================================================
  !     ! sort the nnz(ia ja a) in the order according to 
  !     !   the new permutation It help me on sending 
  !     ======================================================

  tsort=MPI_Wtime()

  !-----------------------------------------------------------------------------
  ! [2.1] Sort the row indices and save the permutation made into "need_Lg"
  !-----------------------------------------------------------------------------
  
  ! sort at an increasing order the row index(iainit) of A and get the permutation order
  ! Warning : IPSORT may fail on really huge matrices. (found by Azzam)
  ! alternative with license issue : KB07AI_mod (HSL)

  istep = 21
  kflag = 2 ! sort in increasing order
  CALL IPSORT(ia,nnza,need_Lg,kflag,iinfo)
  ! CALL KB07AI_mod(ia,nnza,need_Lg)  ! alternative with license issue.
  If (iinfo /= 0 ) iinfo = -iinfo
  If (iinfo <  0 ) Goto 9999
  
  !-----------------------------------------------------------------------------
  ! [2.2] Apply the permutation to the column indices
  !-----------------------------------------------------------------------------
  
  istep = 22
  Call XMPH_iperm(ja, need_Lg, nnza, iinfo)
  If (iinfo <  0 ) Goto 9999

  !-----------------------------------------------------------------------------
  ! [2.3] Apply the permutation to the values
  !-----------------------------------------------------------------------------
  
  istep = 23
  Call XMPH_xperm( a, need_Lg, nnza, iinfo)
  If (iinfo <  0 ) Goto 9999

  !-----------------------------------------------------------------------------
  ! [3] Compute the compressed row vector and the array "need_Lg"
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! [3.1] reset the arrays to compute
  !-----------------------------------------------------------------------------

  Do i=1,nnza
     need_Lg(i)=0
  End Do

  Do i= 1,ndof+1
     iacsr(i) = 0
  End Do

  
  !-----------------------------------------------------------------------------
  ! [3.2] count the number of occurence of each row ("iacsr") and fill "need_Lg"
  !-----------------------------------------------------------------------------

  do i=1,nnza

     ! Count the number of occurence of each row
     
     ivp=ia(i)
     iacsr(ivp+1)= iacsr(ivp+1) + 1


     ! Fill "need_Lg"

     ! -----------------------------------
     ! For Lagrange stuff Azzam:15/04/2009
     ! find non zero diagonal --> the 
     ! remaining zero entries are Lagrange
     ! -----------------------------------

     if(ivp.eq.ja(i)) then ! this is a diagonal entry
        if(a(i).ne.XMPH_FLOATZERO) then ! non zero diagonal
           !C if(abs(a(i)).le.dropLG) then ! non zero diagonal
           need_Lg(ivp)=1
        endif
     endif


  enddo
  
  !-----------------------------------------------------------------------------
  ! [3.3] From the number of occurence of each row, form the compressed vector.
  !-----------------------------------------------------------------------------

  ! add this nb of nnz by row to find the starting point of each 
  iacsr(1)=1
  do i=2,ndof+1
     iacsr(i)= iacsr(i) + iacsr(i-1)
  enddo

  tsort=MPI_Wtime()-tsort

  !-----------------------------------------------------------------------------
  ! [4] End
  !-----------------------------------------------------------------------------

  ! Save data
  sm_global%fmt =  SM_FMT_CSR
  sm_global%csr => iacsr

  ! Print error/warnings
9999 Continue
  If ( iinfo /=  0 ) Then
     
     If ( iinfo > 0) msg_class = MSG_WARNING
     If ( iinfo < 0) msg_class = MSG_ERROR
     
     Select Case(istep) 
     Case(1); Call mph_logWithInfo (msg_class,ndof+1,Trim(rname)//&
          " while allocating the compressed row, size = ")
     Case(21); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
          " while sorting the row indices ")
     Case(22); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
          " while applying permutation to the column indices ")
     Case(23); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
          " while applying permutation to the values ")
     End Select

  End If

  ! Set the return code
  If ( iinfo == 0 ) info =  0
  If ( iinfo <  0 ) info = -istep
  If ( iinfo >  0 ) info = +istep

end Subroutine XMPH_Part_sort_GlobalMatrix_CSR

! [+] routine : XMPH_Part_compute_lagrange ------------------------------------------
!
!> compute informations related to lagrange specific treatment (domLg)
!!
!! @param [in] tree  
!!        The partitioning tree 
!!
!! @param [in] need_Lg  
!!        Integer array specifying if a row/column needs Lagrange treatment.
!!
!! @param [out] domLg  
!!        Integer array of size nbdom + 1.
!!        It specifies the number of lagranges in each domains and on the interface.
!!        - (1..nbdom  ) = number of Lagranges (rows/columns) in interior 1..nbdom.
!!        - (nbdom+1   ) = number of Lagranges (rows/columns) in interfaces.
!!
!! @param [out] info
!!        The return code of the routine.
!!
!! @author Azzam Haidar
!! @author Yohan Lee-tin-yien
Subroutine XMPH_Part_compute_lagrange               &  
     (tree, need_Lg, &
     & domLg, info )
  
  !* Modules *!
  
  use MPH_part_type
  Use mph_error_mod
  Implicit none

  !* Arguments *!

  Type(maphys_binary_tree_t) , intent(in) :: tree
  MPH_INT  , pointer , intent(in  ) :: need_Lg(:)
  Integer, pointer    , intent(out) :: domLg (:)
  Integer             , intent(out) :: info

  !* Local variables *!

  ! Scalars

  Integer    :: nbdom        ! components of tree
  Integer    :: totinterface ! components of tree
  Integer    :: nproc        ! components of tree
  Integer    :: nbsep        ! components of tree
  Integer    :: firstrow     ! the first row/col of a domain
  Integer    :: lastrow      ! the last row/col of a domain
  Integer    :: domid        ! a selected domain,
  Integer    :: procid       ! a selected process
  Integer    :: sepid        ! a selected separator
  MPH_INT :: i            ! a counter

  ! Arrays
  Integer, pointer    :: domptr (:)             ! component of tree 

  ! Derived types
  Type(maphys_binary_node_t), pointer :: azz(:) ! component of tree 

  !- End of header -------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! [1] Init
  !-----------------------------------------------------------------------------

  
  nbsep        =  tree%nbsep
  nbdom        =  tree%nbsep + 1
  nproc        =  nbdom
  totinterface =  tree%totinterface
  azz          => tree%node
  domptr       => tree%domptr

  Allocate(domLg(nproc+1), STAT= info )
  CHCKASSRT( info == 0, info )
  If (info < 0) Return

  Do i=1, nproc+1
     domLg(i)  = 0
  End Do

  !-----------------------------------------------------------------------------
  ! [2] Fill domLg
  !-----------------------------------------------------------------------------

  do domid=1,nproc
     procid = domid-1
     sepid  = domptr(domid)
     if(azz(sepid)%rdomid.eq.procid) then 
        firstrow = azz(sepid)%rgst
        lastrow  = azz(sepid)%rged
        do i=firstrow,lastrow
           if(need_Lg(i).eq.0)domLg(domid)=domLg(domid)+1 
        enddo
     else 
        firstrow = azz(sepid)%lgst
        lastrow  = azz(sepid)%lged
        do i=firstrow,lastrow
           if(need_Lg(i).eq.0)domLg(domid)=domLg(domid)+1 
        enddo
     endif
  enddo
  do i=1,totinterface
     if(need_Lg(i).eq.0)domLg(nproc+1)=domLg(nproc+1)+1
  enddo
     
End Subroutine XMPH_Part_compute_lagrange

#if MAPHYS_DEBUG
! [+] routine : XMPH_Part_compute_lagrange_debug -------------------------------------
!
!> Export lagrange data for checking (with matlab for example)
!!
!! @param[in ] ndof     Order of the global matrix
!! @param[in ] met_perm Array specifying the permutation on row/columns
!! @param[in ] need_Lg  Array specifying Lagrange treatment
!! @param[in ] tree     Partitioning binary tree
!! @param[out] info     Return code
!!
Subroutine XMPH_Part_compute_lagrange_debug(ndof, metperm, need_Lg, tree, info )

  !* Module(s) *!
  Use MPH_part_type

  Use mph_log_mod
  Implicit none


  !* Arguments *!
  integer                    , intent(in) :: ndof
  integer, pointer           , intent(in) :: metperm   (:)
  integer, pointer           , intent(in) :: need_Lg  (:)
  type(maphys_binary_tree_t) , intent(in) :: tree
  integer                    , intent(out) :: info

  !* Local variable *!

  ! scalar
  Integer :: debug
  Integer :: i,j,k, l 
  Integer :: istep, iinfo , msg_class
  Integer :: nbsep                ! component of "tree"
  Integer :: nbdom                ! component of "tree"
  Integer :: totinterface         ! component of "tree"
  Integer :: nproc                ! component of "tree"
  integer :: ios1                 ! file related
  integer :: firstrow, lastrow    ! selected element
  integer :: domid, procid, sepid ! selected element

  ! array
  Integer, allocatable :: intind(:) ! computed & written data
  Integer, allocatable :: lgind(:)  ! computed & written data
  Integer, pointer    :: domptr (:) ! component of "tree"

  ! type
  type(maphys_binary_node_t), pointer :: azz(:)! component of "tree"

  ! string
  Character(len=MAPHYS_STRL), Parameter :: rname = "part_compute_lagrange_debug"
  Character(len=MAPHYS_STRL) :: pbname = "ana_lagrange_debug" ! file related
  Character(len=MAPHYS_STRL) :: filename                      ! file related
  Character(len=MAPHYS_STRL) :: pass                          ! file related

  !- End of header -------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! [1] Init
  !-----------------------------------------------------------------------------

  nbsep        =  tree%nbsep
  nbdom        =  tree%nbsep + 1
  nproc        =  nbdom
  totinterface =  tree%totinterface
  azz          => tree%node
  domptr       => tree%domptr

  !-----------------------------------------------------------------------------
  ! [2] Export separators
  !-----------------------------------------------------------------------------

  debug=0
  if(debug.eq.1)then
     write(unit=filename,FMT='(A,"_separator")') trim(pbname)
     open(unit=18,file=filename,status="REPLACE",iostat=ios1)
     do i=1, nbsep 
        do j=azz(i)%sepst,azz(i)%seped
           write(unit=18,FMT='(I16)')metperm(j)
        enddo
     enddo
     close(unit=18)
  endif

  !-----------------------------------------------------------------------------
  ! [3] Compute & export intind and lgind
  !-----------------------------------------------------------------------------

  istep = 3
  ALLOCATE(intind(ndof),STAT=iinfo)
  If (iinfo > 0) iinfo = -iinfo
  If (iinfo < 0) Goto 9999

  ALLOCATE(lgind(ndof),STAT=iinfo)
  If (iinfo > 0) iinfo = -iinfo
  If (iinfo < 0) Goto 9999

  do domid=1,nproc
     procid = domid-1
     sepid  = domptr(domid)
     k=0
     l=0
     if(azz(sepid)%rdomid.eq.procid) then 
        firstrow = azz(sepid)%rgst
        lastrow  = azz(sepid)%rged
        do i=firstrow,lastrow
           if(need_Lg(i).eq.0)then
              k=k+1
              lgind(k)=i
           else
              l=l+1
              intind(l)=i
           endif
        enddo
     else 
        firstrow = azz(sepid)%lgst
        lastrow  = azz(sepid)%lged
        do i=firstrow,lastrow
           if(need_Lg(i).eq.0)then
              k=k+1
              lgind(k)=i
           else
              l=l+1
              intind(l)=i
           endif
        enddo
     endif


     if(debug.eq.1)then
        !        write interior indices to a file for matlab
        write(unit=pass,FMT='(I5)')domid
        write(unit=filename,FMT='(A,"_intind_",A)')&
             trim(pbname),adjustl(pass)
        open(unit=18,file=filename,status="REPLACE",iostat=ios1)
        do j=1,l 
           write(unit=18,FMT='(I16)')intind(j)
        enddo
        close(unit=18)
        !        write lagrange indices to a file for matlab
        write(unit=pass,FMT='(I5)')domid
        write(unit=filename,FMT='(A,"_lgind_",A)')&
             trim(pbname),adjustl(pass)
        open(unit=18,file=filename,status="REPLACE",iostat=ios1)
        do j=1,k 
           write(unit=18,FMT='(I16)')lgind(j)
        enddo
        close(unit=18)
     endif
  enddo

  !-----------------------------------------------------------------------------
  ! [4] End
  !-----------------------------------------------------------------------------

  ! Print error/warnings

9999 Continue
  If ( iinfo /=  0 ) Then
     
     If ( iinfo > 0) msg_class = MSG_WARNING
     If ( iinfo < 0) msg_class = MSG_ERROR
     
     Select Case(istep) 
     Case(1); Call mph_logWithInfo (msg_class,nproc+1,Trim(rname)//&
          " while allocating the array need_Lg, size = ")
     End Select
  End If

  ! Set the return code

  If ( iinfo == 0 ) info =  0
  If ( iinfo <  0 ) info = -istep
  If ( iinfo >  0 ) info = +istep

  ! Free memory
  If (Allocated(lgind )) Deallocate(lgind)
  If (Allocated(intind)) Deallocate(intind)

end Subroutine XMPH_Part_compute_lagrange_debug
#endif


! [+] routine : XMPH_Part_compute_domains_nnz ---------------------------------------
!
!> compute all number of non-zero entry for each domain
!! ie: nnzallinterface, domintnnz, domintrfnnz 
!! @param [in] tree  
!!        The partitioning tree 
!!
!! @param [in] sm_global
!!        The matrix (previously sorted in CSR)
!!
!! @param [in,out] domains  
!!        Specifies the domain decomposition. 
!!        On output, append the computation of those components :
!!        - nnzallinterface
!!        - domintnnz
!!        - domintrfnnz
!!        .
!! @param [out] info
!!        The return code of the routine.
!!
!! @author Azzam Haidar
!! @author Yohan Lee-tin-yien
!!
Subroutine XMPH_Part_compute_domains_nnz( tree, sm_global, domains, info )

  !* Modules *!
  Use MPH_part_type
  Use XMPH_sparse_matrix_mod
  Use mph_log_mod
#if MAPHYS_DEBUG
  Use mph_dbg_mod
#endif

  Implicit none

  !* Arguments *!
  Type(maphys_binary_tree_t) , Intent(in) :: tree
  Type(XMPH_sparse_matrix_t) , Intent(in) :: sm_global
  Type(maphys_domains_t), Intent(inout)   :: domains
  Integer                                 :: info

  !* Local variables *!

  ! Scalar
  Integer    :: istep           ! internal step
  Integer    :: iinfo           ! called routine's status
  Integer    :: msg_class       ! class of message
  Integer    :: firstrow        ! the first row/col of a domain
  Integer    :: lastrow         ! the last row/col of a domain
  Integer    :: domid           ! a selected domain,
  Integer    :: procid          ! a selected process
  Integer    :: sepid           ! a selected separator
  Integer    :: nbdom           ! component of "tree"
  Integer    :: nproc           ! component of "tree"
  MPH_INT :: totinterface    ! component of "tree"
  MPH_INT :: nnzallinterface ! component of "domains"

  ! Array
  Integer    , Pointer :: domptr (:)     ! component of "tree"
  MPH_INT , Pointer :: domintnnz(:)   ! component of "domains"
  MPH_INT , Pointer :: domintrfnnz(:) ! component of "domains"
  MPH_INT , Pointer :: iacsr(:)       ! component of "sm_global"

  ! String
  Character(len=MAPHYS_STRL), Parameter :: rname = "XMPH_Part_Compute_domains_nnz"

  ! Derived type
  Type(maphys_binary_node_t), pointer :: azz(:) ! component of "tree"


  ! -- selected elements

  ! End of header --------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! [1] Init
  !-----------------------------------------------------------------------------

  !
  istep = 11
  iacsr => sm_global%csr

  nbdom        =  tree%nbsep + 1
  nproc        = nbdom
  totinterface =  tree%totinterface
  azz          => tree%node
  domptr       => tree%domptr

  !
  istep = 12
  Allocate(domintnnz(nproc)  , STAT=iinfo)
  If (iinfo > 0) iinfo = -iinfo
  If (iinfo <  0) Goto 9999

  Allocate(domintrfnnz(nproc), STAT=iinfo)
  If (iinfo > 0) iinfo = -iinfo
  If (iinfo < 0) Goto 9999

  Do domid=1,nproc
     domintnnz(domid)  = 0
     domintrfnnz(domid)  = 0
  End do


  !-----------------------------------------------------------------------------
  ! [2] Compute the domain nnz.
  !-----------------------------------------------------------------------------
  !
  ! count the interior_row nnz by domain
  ! here interior row are formed from 
  ! [all interior + interior interface connexion]
  !

  do domid=1,nproc
     procid = domid-1
     sepid  = domptr(domid)
     if(azz(sepid)%rdomid.eq.procid) then 
        firstrow = azz(sepid)%rgst
        lastrow  = azz(sepid)%rged
        domintnnz(domid) = iacsr(lastrow+1)-iacsr(firstrow)
        if(domid.eq.1)nnzallinterface=iacsr(firstrow)-1
     else 
        firstrow = azz(sepid)%lgst
        lastrow  = azz(sepid)%lged
        domintnnz(domid) = iacsr(lastrow+1)-iacsr(firstrow)
        if(domid.eq.1)nnzallinterface=iacsr(firstrow)-1
     endif
  enddo

  !-----------------------------------------------------------------------------
  ! [3] End
  !-----------------------------------------------------------------------------

  ! Save data
  domains%nnzallinterface = nnzallinterface
  domains%domintrfnnz     => domintrfnnz
  domains%domintnnz       => domintnnz

#if MAPHYS_DEBUG
  ! Export data
  Call MPH_dbg_init
  Call MPH_dbg_set_file("domains_nnz.txt")
  Write(dbg_unit,*) "nnzallinterface =", domains%nnzallinterface
  Write(dbg_unit,*) "domintrfnnz     =", domains%domintrfnnz (1:nproc)
  Write(dbg_unit,*) "domintnnz       =", domains%domintnnz   (1:nproc)
  Call MPH_dbg_exit
#endif

  ! Print error/warnings
9999 Continue
  If ( iinfo /=  0 ) Then
     
     If ( iinfo > 0) msg_class = MSG_WARNING
     If ( iinfo < 0) msg_class = MSG_ERROR
     
     Select Case(istep) 
     Case(12); Call mph_logWithInfo (msg_class,nproc,Trim(rname)//&
          " while allocating domintnnz or domintrfnnz, size = ")
     End Select

  End If

  ! Set the return code
  If ( iinfo == 0 ) info =  0
  If ( iinfo <  0 ) info = -istep
  If ( iinfo >  0 ) info = +istep


End Subroutine XMPH_Part_compute_domains_nnz


End Module XMPH_part_builddomains_mod
