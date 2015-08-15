! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"
! [+] module : XMPH_part_distmatrix_mod --------------------------------
!
!> module to distribute the matrix
!!
Module XMPH_part_distmatrix_mod

  !* Module(s) *!
  Use mph_error_mod
  
  !* No implicit typing *!
  Implicit None

  !* Private constants *!
  Character(len=MAPHYS_STRL), Private, Parameter :: &
       FLNAME= "XMPH_ARITHmph_part_distmatrix_mod.F90"

  !* Defined routine(s) *!
  Public :: XMPH_PART_DistGlobalMatrix
  Public :: XMPH_PART_compute_interface_weight

Contains

  ! [+] routine : XMPH_PART_DistGlobalMatrix ---------------------------------
  !
  !> Distribute the Global matrix.
  !!
  !! which means, send to each processor their local system + domain
  !! informations and on processor 0, save paritionning informations to
  !! split the future second member.
  !! 
  !!-----
  !!
  !! @param[in     ] me            the MPI rank 
  !! @param[in     ] comm          the MPI communicator
  !! @param[in     ] all_domains   the datas about domains partitioning
  !!
  !! @param[in,out ] sm_global     the matrix to distribute
  !! @param[in,out ] graph_global  the global matrix graph
  !! @param[in,out ] ana_timing    timings of the analysis phase
  !!
  !! @param[out    ] part_rhs  the partioning informations 
  !!                               neccessary to split the second member
  !! @param[out    ] one_domain    the local domain description
  !! @param[out    ] sm_local      the local domain matrix
  !! @param[out    ] info          the subroutine status
  !! 
  !!-----
  !!
  !! @warning
  !! the subroutine deallocate the memory of its inputs;
  !! hence the intent(inout) flags for INPUTS.
  !! 
  !! @todo finish reformating comments
  !! @todo write algorithm comments
  !! @todo move [dbg.1.4] in another routine
  !!
  !! @author Azzam Haidar
  !! @author Yohan Lee-tin-yien
  !!
  !! @version 0.2b Yohan : subroutine created & renamed
  !! @version 0.2a Azzam : update of read_partmatices from Azzam (new comments)
  !! @version 0.1  Azzam : Azzam's initial "read_partmatrices"  
  !!
  Subroutine  XMPH_PART_DistGlobalMatrix  &
       (me, comm, all_domains,              &
       sm_global, graph_global, ana_timing, &
       part_rhs, one_domain, sm_local,      &
       info )

    !* Module(s) *!
    
    Use XMPH_sparse_matrix_mod
    Use MPH_part_type
    Use mph_log_mod
    Use mph_error_mod
    implicit none
    include 'mpif.h'

    !* Arguments *!

    integer                     , intent(in)    :: me      ! MPI rank
    integer                     , intent(inout) :: comm    ! MPI communicator
    type(maphys_domains_t)      , intent(inout) :: all_domains

    type(XMPH_sparse_matrix_t) , intent(inout) :: sm_global
    type(maphys_matrix_graph_t) , intent(inout) :: graph_global
    real(kind=8), intent(inout) :: ana_timing(ANA_TIMING_SIZE)

    type(maphys_rhs_partition_t) , intent(out) :: part_rhs
    type(maphys_domain_t   )    , intent(out) :: one_domain
    type(XMPH_sparse_matrix_t) , intent(out) :: sm_local
    integer                     , intent(out) :: info

    !* Local variables *!

    ! contants
    Integer            :: master  = 0
    integer, parameter :: RHSWAY  = 2
    integer, parameter :: WRONG_INT_VALUE  = -1
    real   , parameter :: WRONG_REAL_VALUE = -1.d0

    ! sm_global association
    MPH_INT             :: nnza
    MPH_INT   , pointer :: iacsr     (:)
    MPH_INT   , pointer :: ia        (:)
    MPH_INT   , pointer :: ja        (:)
    XMPH_FLOAT , pointer :: a         (:)

    ! all_domains association
    integer :: combivtxpcsz
    integer :: nproc
    integer :: maxprocintrf
    MPH_INT       :: nnzallinterface

    integer, pointer :: domstptr    (:)

    integer, pointer :: domintnnz     (:)
    integer, pointer :: domintrfnnz   (:)
    integer, pointer :: domintdof     (:)
    integer, pointer :: domintrfdof   (:)

    integer, pointer  :: intrfindices (:)
    integer, pointer  :: intrfproc    (:)

    integer, pointer :: vtxpos        (:)
    integer, pointer :: vtxweight     (:)

    integer, pointer :: domLg         (:)

    integer :: totinterface
    integer, pointer :: metperm           (:)

    ! sm_local association
    integer               :: symtype
    integer               :: myndof
    MPH_INT            :: mynnza
    MPH_INT  , pointer :: myiacsr  (:)
    MPH_INT  , pointer :: myia     (:)
    MPH_INT  , pointer :: myja     (:)
    XMPH_FLOAT, pointer :: mya      (:)

    ! graph_global association
    integer :: ndof
    integer :: maxadj
    integer, pointer :: adjncy        (:)


    ! paritioning association
    ! integer, pointer :: metperm           (:)

    ! integer ::  combivtxpcsz
    integer, pointer :: procintdisp       (:)
    integer, pointer :: procintrfdisp     (:)
    ! integer, pointer :: domLg             (:)
    ! integer, pointer :: domintdof         (:)
    ! integer, pointer :: domintrfdof       (:)
    integer, pointer :: domlogicintrfdof  (:)
    integer, pointer :: scatindices       (:)
    integer, pointer :: scatlogicindices  (:)

    ! one_domain association
    integer :: mynbvi
    integer :: myndofinterior
    integer :: myndofintrf
    integer :: myndoflogicintrf
    integer :: mysizeIntrf
    integer :: myint_Lg
    integer :: gballndof
    integer :: gballintrf
    integer :: lenindintrf

    integer, pointer :: myindexVi      (:)
    integer, pointer :: myptrindexVi   (:)
    integer, pointer :: myindexintrf   (:)
    integer, pointer :: myinterface    (:)
    integer, pointer :: mylogicintrf   (:)
    integer, pointer :: gbtoloc        (:)

    !-- internal status (errors & checks)
    Integer :: infompi  ! status for calls to MPI routines
    Integer :: iinfo    ! status for calls to routines
    Integer :: ierr
    Integer :: ierr1
    Integer :: istep     ! internal step in the routine
    Integer :: msg_class ! class of a message.
    Character(len=MAPHYS_STRL) :: msg ! message

    !while checking duplicate entries
    integer :: thisindex
    integer :: thisproc

    integer :: previndex
    integer :: prevproc

    !-- timings
    real(kind=8) :: tmetis
    real(kind=8) :: tsort
    real(kind=8) :: ttmp

    real(kind=8) :: tdata   ! time to construct & send domain structure
    real(kind=8) :: tscat   ! time to scatter the global matrix
    real(kind=8) :: tloc    ! time to convert to local indices
    real(kind=8) :: tparti  ! time to partition matrix

    !-- temporary data
    integer :: i     ! counter on processus
    integer :: j,k,l ! counters
    integer :: t,u,w ! counters

    integer :: pt    ! counter

    integer :: m     ! counter on domains
    integer :: intcntproc !azz: for each neighbor proc, this is the number of sharing interface

    integer :: st    ! starting index
    integer :: ed    ! ending   index
    integer :: domst ! domain starts at this index
    integer :: domed ! domain ends   at this index

    integer :: pcnb    ! processus number
    integer :: procid  ! a processus (processus indentifier)
    integer :: domid   ! a domain    (domain    indentifier)
    integer :: scatpt
    integer :: scatptlogic

    integer :: nbVi   ! number of neighbors (French : "nombre de Voisins" )
    integer :: vtxid  ! selected vertex

    integer :: rowsiz ! size of a local matrix row
    ! used for calls to isort
    integer :: kflag
    integer :: nothing

    ! only useful while computing the local interface structure
    integer          :: intrfsize
    integer, pointer :: listinterface  (:)
    integer, pointer :: indexIntrf     (:)
    integer, pointer :: indexIntrfproc (:)
    integer, pointer :: indexVi        (:)
    integer, pointer :: ptrindexVi     (:)

    integer :: right, left
    integer :: rint,  lint

    ! only useful while computing the local interfae values
    integer, pointer :: belongingto    (:)
    integer, pointer :: tempsepi       (:)
    integer, pointer :: tempsepj       (:)
    XMPH_FLOAT, pointer :: tempsepval(:)
    integer :: jcol ! selected column
    integer :: jptr ! pointer on columns

    ! used while sending domain structure
    MPH_INT :: nnzaintrf      ! number of non-zeros on this interface
    MPH_INT :: nnzainterior   ! number of non-zeros on this interior
    MPH_INT :: allsendnnz     !

    integer, parameter :: GINFOSIZE = 20 ! is 15 not enougth ?
    MPH_INT :: infoproc(GINFOSIZE)   ! buffer storing local system sizes (send)
    integer    :: tagg           ! starting tag index for MPI communications

    ! used while recveiving domain structure
    INTEGER    :: stat(MPI_STATUS_SIZE) ! status in calls to MPI routines
    MPH_INT :: myinfoproc(GINFOSIZE) ! buffer storing local system sizes (recv)
    MPH_INT :: mynnzaintrf      ! number of non-zeros on this interface
    MPH_INT :: mynnzainterior   ! number of non-zeros on this interior

    ! used while sending interior rows
    MPH_INT, pointer :: procintcnt (:)
    MPH_INT, pointer :: procnzdisp (:)

    ! used while switching to local indexing
    MPH_INT :: jvp

    ! used for debuging purpose
    integer :: debug 
    integer :: p, nbsep

    integer :: unitnb5 = 31
    integer :: ios1
    character(MAPHYS_STRL) :: pbname = "PART_DistGlobalMatrix_dbg"
    character(MAPHYS_STRL) :: filename
    character(MAPHYS_STRL) :: pass


    !-- others
    MPH_INT :: BUGMUMP ! ? = mysizeIntrf + myint_Lg ?

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [0] Initialize local variables
    !---------------------------------------------------------------------------
    ! 
    ! Init errors
    infompi = 0
    iinfo   = 0
    ierr    = 0
    ierr1   = 0

    ! Init counters
    allsendnnz = 0

    ! Nullify internal allocated memory

    ! allocated in section [2.1]
    Nullify(procintcnt )
    Nullify(procintdisp)
    Nullify(procintrfdisp)
    Nullify(procnzdisp )
    Nullify(indexIntrf)
    Nullify(listinterface)

    ! Init outputs
    Call XMPH_sm_nullify(sm_local,iinfo)

    part_rhs%rhsway = WRONG_INT_VALUE 
    part_rhs%combivtxpcsz = WRONG_INT_VALUE                
    Nullify(part_rhs%procintdisp      )      
    Nullify(part_rhs%procintrfdisp    )   
    Nullify(part_rhs%domLg            )            
    Nullify(part_rhs%domintdof        )        
    Nullify(part_rhs%domintrfdof      )     
    Nullify(part_rhs%metperm          )          
    Nullify(part_rhs%scatindices      )
    Nullify(part_rhs%scatlogicindices )
    Nullify(part_rhs%domlogicintrfdof )

    combivtxpcsz = WRONG_INT_VALUE                
    Nullify(procintdisp      )      
    Nullify(procintrfdisp    )   
    Nullify(domLg            )            
    Nullify(domintdof        )        
    Nullify(domintrfdof      )     
    Nullify(metperm          )          
    Nullify(scatindices      )
    Nullify(scatlogicindices )
    Nullify(domlogicintrfdof )

    ! Associate local variables to in and inouts variables
    if ( me == 0) then ! master have the data
       !-- sm_global
       symtype      =  sm_global%sym
       nnza         =  sm_global%nnz
       ia           => sm_global%i
       ja           => sm_global%j
       a            => sm_global%v
       iacsr        => sm_global%csr
       !-- all_domains
       combivtxpcsz =  all_domains%combivtxpcsz
       nproc        =  all_domains%nbdom
       maxprocintrf =  all_domains%maxprocintrf
       nnzallinterface = all_domains%nnzallinterface

       domintdof    => all_domains%domintdof
       domintrfdof  => all_domains%domintrfdof

       domintnnz    => all_domains%domintnnz
       domintrfnnz  => all_domains%domintrfnnz

       intrfindices => all_domains%intrfindices
       intrfproc    => all_domains%intrfproc

       vtxweight => all_domains%vtxweight
       vtxpos    => all_domains%vtxpos

       domLg     => all_domains%domLg
       domstptr  => all_domains%domstptr

       totinterface =  all_domains%totinterface
       metperm      => all_domains%metperm

       !-- graph_global
       ndof         =  graph_global%ndof
       maxadj       =  graph_global%maxadj
       adjncy       => graph_global%adjncy

       ! ana_timing
       tparti  = ana_timing(7)

    else ! Give default values to others 

       ! sm_global
       symtype      = WRONG_INT_VALUE
       nnza         = WRONG_INT_VALUE
       nullify(ia)
       nullify(ja)
       nullify(a)
       nullify(iacsr)

       ! all_domains
       combivtxpcsz = WRONG_INT_VALUE
       nproc        = WRONG_INT_VALUE
       maxprocintrf = WRONG_INT_VALUE
       nnzallinterface = WRONG_INT_VALUE
       totinterface = WRONG_INT_VALUE

       nullify(domintdof)
       nullify(domintrfdof)

       nullify(domintnnz)
       nullify(domintrfnnz)

       nullify(intrfindices)
       nullify(intrfproc)

       nullify(vtxweight)
       nullify(vtxpos)

       nullify(domLg)
       nullify(domstptr)

       nullify(metperm)

       ! graph_global
       ndof         = WRONG_INT_VALUE
       maxadj       = WRONG_INT_VALUE
       nullify(adjncy)

       ! ana_timing
       tparti  = WRONG_REAL_VALUE

    endif


    !---------------------------------------------------------------------------
    ! [1] Compute & Send or Receive the domain's structure 
    !---------------------------------------------------------------------------

    if ( me == 0) then
       tdata = MPI_Wtime()

       !------------------------------------------------------------------------
       ! [1.0.0] Initialize & allocate temporary data
       !------------------------------------------------------------------------

       !> @par Description
       !! on Master proc 0 :
       !! for each proc construct the data structure and send it.
       !! this part can be done in parallel by each processor so
       !! the master have to BCAST [intrfindices intrfproc] to all
       !! the processor and then each one compute its data structure
       !!
       !! @Attention  
       !! [intrfindices intrfproc] are sorted according
       !! to the indices. intrfindices is sorted and carry intrfproc along
       !!
       !!

       !azzam>
       !C      IF(USINGBSEND.EQ.0)THEN
       !C*        Attach a buffer big enough for sending the indices, in bytes
       !C         buffsiz=2*(combivtxpcsz+nproc)
       !C         ALLOCATE(combuff(buffsiz + 10*MPI_BSEND_OVERHEAD))
       !C         sizebuff = (buffsiz + 10*MPI_BSEND_OVERHEAD)*MPH_INTBYTESIZE
       !C         CALL MPI_BUFFER_ATTACH(combuff, sizebuff, infompi)
       !C      ENDIF
       !azzam<

       !> @par Strategy description :
       !!
       !! - way1 to generate the scatterv indices: 
       !!   for each proc send the corresponding interface RHS
       !!   and interior RHS
       !!
       !! - way2 to generate the scatterv indices: 
       !!   for each proc send the corresponding logical interface RHS 
       !!   and interior RHS
       istep = 100
       !way1
       IF(RHSWAY.EQ.1)Then
          Allocate(scatindices(combivtxpcsz), STAT= iinfo)
          If (iinfo > 0) iinfo = -iinfo 
          If (iinfo < 0) Goto 9999
       End IF

       !way2
       IF(RHSWAY.EQ.2)Then

          ! scatlogicindices
          ALLOCATE(scatlogicindices(totinterface), stat=iinfo)
          If (iinfo > 0) iinfo = -iinfo 
          If (iinfo < 0) Goto 9999

          Do i=1,totinterface
             scatlogicindices(i)=  0
          End Do

          ! domlogicintrfdof
          ALLOCATE(domlogicintrfdof(nproc), stat=iinfo)
          If (iinfo > 0) iinfo = -iinfo 
          If (iinfo < 0) Goto 9999
          Do i=1,nproc
             domlogicintrfdof(i)=  0
          End Do

       End IF

       ALLOCATE(listinterface(totinterface), STAT=iinfo)
       If (iinfo > 0) iinfo = -iinfo 
       If (iinfo < 0) Goto 9999

       ALLOCATE(indexIntrf(combivtxpcsz), STAT=iinfo)
       If (iinfo > 0) iinfo = -iinfo 
       If (iinfo < 0) Goto 9999

       ALLOCATE(indexIntrfproc(combivtxpcsz), STAT=iinfo)
       If (iinfo > 0) iinfo = -iinfo 
       If (iinfo < 0) Goto 9999

       ALLOCATE(indexVi(nproc-1), STAT=iinfo)
       If (iinfo > 0) iinfo = -iinfo 
       If (iinfo < 0) Goto 9999

       ALLOCATE(ptrindexVi(nproc), STAT=iinfo)
       If (iinfo > 0) iinfo = -iinfo 
       If (iinfo < 0) Goto 9999

       ALLOCATE(tempsepi(maxprocintrf+1), STAT=iinfo)
       If (iinfo > 0) iinfo = -iinfo 
       If (iinfo < 0) Goto 9999

       ALLOCATE(tempsepj(maxprocintrf*(maxadj+1)), STAT=iinfo)
       If (iinfo > 0) iinfo = -iinfo 
       If (iinfo < 0) Goto 9999

       ALLOCATE(tempsepval(maxprocintrf*(maxadj+1)), STAT=iinfo)
       If (iinfo > 0) iinfo = -iinfo 
       If (iinfo < 0) Goto 9999

       ALLOCATE(belongingto(ndof), STAT=iinfo)
       If (iinfo > 0) iinfo = -iinfo 
       If (iinfo < 0) Goto 9999


       !way1 or way2
       adjncy(:)   = 0 ! to notify if a entries A(i,j) is already treated by a processor or not.
       ! I need just the first nnza value of adjncy.
       ! SEE BELOW generate (ia ja a) part for details
       scatpt=0
       scatptlogic=0

       do i=1,nproc
          intrfsize    = 0
          pt           = 0
          pcnb         = i-1
          indexIntrf(:)    = -1
          indexIntrfproc(:)= -1
          listinterface(:) = -1
          tempsepi(:)   = 0
          tempsepj(:)   = 0
          tempsepval(:) = XMPH_FLOATZERO

          !---------------------------------------------------------------------
          ! [1.1] Compute the domain's structure 
          !---------------------------------------------------------------------

          !azz>
          ! ATTENTION [intrfindices intrfproc] are sorted according
          ! to intrfindices. So the intrfindices is sorted and carry intrfproc along
          !
          ! first find all the interface indices of this proc
          ! and their corresponding neighbor processors
          ! listinterface : will contain list of interface of proc pcnb (not repeted indices)
          ! indexIntrf    : will contain the indices of the interface repeted
          ! intrfproc     : will contain the corresponding neighbor
          !                   procId of the indices of indexIntrf
          ! scatindices or       : will contain the indices of the interface of all proc and sorted by proc.
          ! scalogictindices     [intrftoproc0 intrftoproc1 intrftoproc2 ...]
          !
          !                     PAY ATTENTION THAT FOR EACH PROC:
          !                      - indices should be as the same order as numbered by the proc itself.
          !                      - indices should stored after permutation. it is easy to scatterv
          !azz<
          !

          !-------------------------------------------------------------
          ! [1.1.1] Check duplicated entries
          !-------------------------------------------------------------

          !
          !   test that there are not a same
          !   index repeated for the same processor.
          !   Such situation appear if matrix entries
          !   are duplicated
          !

          previndex = intrfindices(1)
          prevproc = intrfproc(1)
          do j=2,combivtxpcsz
             thisindex = intrfindices(j)
             thisproc  = intrfproc(j)
             if(thisindex.eq.previndex)then
                if(thisproc.eq.prevproc)then
                   istep = 111
                   iinfo = -1
                   Goto 9999
                else
                   prevproc = thisproc
                endif
             else
                previndex = thisindex
                prevproc  = thisproc
             endif
          enddo

          !-------------------------------------------------------------
          ! [1.1.2] Compute indexIntrf indexIntrfproc + counters
          !-------------------------------------------------------------

          do j=1,combivtxpcsz
             if(intrfproc(j).eq.pcnb) then !its my interface so search neighbor proc of this index

                !---------------------------------------------
                ! [1.1.2.0] Init variables & increment counters
                !---------------------------------------------

                thisindex=intrfindices(j)
                intrfsize=intrfsize+1
                listinterface(intrfsize)=thisindex

                !---------------------------------------------
                ! [1.1.2.1] way1 to generate the scatterv indices
                !---------------------------------------------

                IF(RHSWAY.EQ.1)THEN
                   scatpt=scatpt+1
                   scatindices(scatpt)=thisindex  ! pay attention to thisindex value
                   ! if it is after or before permutation.
                   ! I need it after to directly scatterv RHS
                   ! otherwise I have to correct when sending and receiving RHS SOL
                ENDIF

                !---------------------------------------------
                ! [1.1.2.2] Treat the right direction
                !---------------------------------------------

                ! right direction: search same indices as thisindex
                ! in the right direction and find each neighbor proc of this index

                do right=j+1,combivtxpcsz
                   rint = intrfindices(right)
                   if(rint.eq.thisindex) then
                      pt=pt+1
                      indexIntrf(pt)=thisindex
                      indexIntrfproc(pt)=intrfproc(right)
                   else
                      goto 323
                   endif
                enddo

                !---------------------------------------------
                ! [1.1.2.3] Treat the left direction
                !---------------------------------------------

                ! left direction: search same indices as thisindex
                ! in the left direction and find each neighbor proc of this index

323             do left=1,j-1
                   lint = intrfindices(j-left)
                   if(lint.eq.thisindex) then
                      pt=pt+1
                      indexIntrf(pt)=thisindex
                      indexIntrfproc(pt)=intrfproc(j-left)
                   else
                      goto 324
                   endif
                enddo
324             continue

             endif ! endif searching neighbor for this index
          enddo ! end loop for all interface indices for proc pcnb

          !---------------------------------------------------
          ! [1.1.2.4] Check data
          !---------------------------------------------------

          ! no interface exist subdomain isolated
          If(pt.eq.0) Goto 973
          ! if(pt.gt.combivtxpcsz) goto 972

          !-------------------------------------------------------------
          ! [1.1.3] Compute the interface
          !-------------------------------------------------------------

          !
          !   now indexIntrf[1:pt] contain the interface indices repeted
          !   and for each value the corresponding neighbor procid
          !   stored in indexIntrfproc[1:pt].
          !
          !   Sort the repeated interface indices of this proc (pcnb)
          !   according to the neighbor processor, that is for each
          !   neighbor proc I will have the corresponding indices
          !   successively. then for each neighbor sort the indices to
          !   keep same order of listing interface between two neighbor
          !   subdomains. then construct the data structure and send it
          !   to proc pcnb.
          !

          !---------------------------------------------------
          ! [1.1.3.1] Sort at an decreasing order 
          !           indexIntrfproc (+ indexIntrf)
          !---------------------------------------------------

          ! Sort at an decreasing order the indexIntrfproc 
          ! sorting the interface indices in the same order

          kflag=-2
          CALL ISORT(indexIntrfproc,indexIntrf,pt,kflag)

          !---------------------------------------------------
          ! [1.1.3.2] Compute indexVi and ptrindexVi
          !---------------------------------------------------

          m=1
          indexVi(m)    = indexIntrfproc(1)
          ptrindexVi(m) = 1
          intcntproc    = 0  ! for each neighbor proc, this is the number of sharing interface

          do l=1,pt
             procid=indexIntrfproc(l)
             if(procid.ne.indexVi(m)) then !new neighbor
                m             = m+1
                indexVi(m)    = procid
                ptrindexVi(m) = ptrindexVi(m-1)+intcntproc

                ! sort the my indexIntrf sharing with proc
                ! [indexVi(m-1)] to force same order of interface
                ! between me and this neighboring processor

                st    = l-intcntproc ! find from where I shall start sorting
                if(st.ne.ptrindexVi(m-1)) goto 974
                kflag = 1
                CALL ISORT(indexIntrf(st),nothing,intcntproc,kflag)
                intcntproc=1
             else
                intcntproc=intcntproc+1
             endif

             if(l.eq.pt) then !last neighbor
                m             = m+1
                ptrindexVi(m) = ptrindexVi(m-1)+intcntproc
                st    = l-intcntproc+1 ! (+1) coz the end find from
                ! where I shall start sorting
                if(st.ne.ptrindexVi(m-1)) goto 974
                kflag = 1
                CALL ISORT(indexIntrf(st),nothing,intcntproc,kflag)
             endif
          enddo
          if(pt.ne.ptrindexVi(m)-1) goto 975
          if(domintrfdof(pcnb+1).ne.intrfsize) goto 973

          !---------------------------------------------------
          ! [1.1.3.3] Generate (ia ja a) of the separator_row
          !       to be sended
          !---------------------------------------------------

          ! put 1 on the corresponding indices of interface and on the
          ! interior indices of this domain then I will search the ja of
          ! the separotor row that correspond to this domain

          ! for each row of the separator_row find the corresponding 
          ! (ja a)
          ! PAY ATTENTION that tempsepi use the same order of indices
          ! as listinterface

          ! For each separator_row:
          !   1-if I am the owner of this vtxid, so I am the last proc
          !     that own this row=vtxid so I should take all the column
          !     that are in my interface even if the indices of one or
          !     more of these columns belong to another proc.
          !     BECAUSE HERE I am the last owner of this row so if I
          !     leave some columns noone can take them.
          !     an another way column are for me when I have to take
          !     their corresponding row. row and col are independant.
          !
          !   2-if I am not the owner of this row, so a processor lowest
          !     than me will take this row with all its interface. So still
          !     remaining column which are not on its interface that is
          !     column which I am their owner and column which I should
          !     be their owner. for this later report to the BUG description.
          !     So I have to take this column in my part.
          !    BUG IN STEP 2
          !

          !   BUG FOR THIS STRATEGY
          !   corrected by the use of adjncy and by giving the lowest proc logically instead of the biggest

          !
          !   For that STEP 2 become:
          !
          !   2-if I am not the owner of this row, so a processor lowest
          !     than me will take this row with all its interface. So still
          !     remaining column which are not on its interface. So columns
          !     that are on my interface and not treated by the owner of
          !     this row and also not treated already by some processor
          !     before me. That is I am the first or the lowest processor
          !     that take this column of this row into account.
          !     So I am the first processor that take this entries [A(106,105)]
          !     into account. For that I use the adjncy vector to indicate
          !     if a column of a row is already treated or not.
          !     NOTE THAT I CAN USE ja INSTEAD OF adjncy AND PUT
          !     ja(w)=-1 BUT I PREFER TO PRESERVE ja

          belongingto(:)=0
          do t=1,intrfsize
             belongingto(listinterface(t)) = 1
          enddo
          domst = domstptr(pcnb+1)
          domed = domstptr(pcnb+1) + domintdof(pcnb+1) -1
          do t=domst,domed
             belongingto(t) = 1
          enddo

          ! put 2 on the interface indices that should belong to a processor lowest than me
          nbvi=m-1
          do w=1,nbvi
             if(indexVi(w).lt.pcnb)then
                do u=ptrindexVi(w),ptrindexVi(w+1)-1
                   belongingto(indexIntrf(u))=2
                enddo
             endif
          enddo
          if(pcnb.eq.(nproc-1))then !last processor
             do t=1,intrfsize
                if(belongingto(listinterface(t)).eq.1)Then
                     Write(msg,'(A,3I10)') 'ERROR LAST PROC HAS LOGIC NTERFACE',&
                     t,intrfsize,listinterface(t)
                     Call mph_log(MSG_ERROR,msg)
                 End if
             enddo
          endif
          !
          jptr        = 0
          tempsepi(1) = 1
          do t=1,intrfsize
             vtxid = listinterface(t)

             !--------------------------------------
             ! [1.1.3.3.1] Save data for scatterv 
             !--------------------------------------

             ! way2 to generate the scatterv indices

             IF(RHSWAY.EQ.2)THEN
                if(belongingto(vtxid).eq.1) then
                   domlogicintrfdof(pcnb+1) = domlogicintrfdof(pcnb+1)+1
                   scatptlogic = scatptlogic +1
                   scatlogicindices(scatptlogic)=vtxid ! pay attention
                   ! to vtxid value if it is after or before
                   ! permutation.  I need it after to directly scatterv
                   ! RHS otherwise I have to correct when sending and
                   ! receiving RHS SOL
                endif
             ENDIF

             if(belongingto(vtxid).eq.1) then  

                !-----------------------------------
                ! [1.1.3.3.1] Handle case "Data belongs to me"
                !-----------------------------------

                ! I am the owner of this row so take all column of
                ! myinterface

                do w=iacsr(vtxid),iacsr(vtxid+1)-1
                   jcol=ja(w)
                   if(belongingto(jcol).ne.0) then     ! interface column
                      jptr  = jptr+1
                      tempsepj(jptr)   = ja(w)
                      tempsepval(jptr) =  a(w)
                      adjncy(w) = -1 ! entry A(vtxid,jcol) is treated
                   endif
                enddo

             elseif(belongingto(vtxid).eq.2) then 

                !-----------------------------------
                ! [1.1.3.3.2] Handle case "Data belongs to my child"
                !-----------------------------------

                ! a proc lowest than me has owned this row.  so take the
                ! remaining column that belong to my interface.

                do w=iacsr(vtxid),iacsr(vtxid+1)-1
                   jcol=ja(w)
                   if(adjncy(w).ne.-1)then             ! this entries A(row col) was not treated. if col is in my interface so I will treat it
                      if(belongingto(jcol).ne.0) then  ! this column is in my interface
                         jptr  = jptr+1
                         tempsepj(jptr)   = ja(w)
                         tempsepval(jptr) =  a(w)
                         adjncy(w) = -1 ! entry A(vtxid,jcol) is treated
                      endif
                   endif
                enddo
             else ! IMPOSSIBLE TO ENTER HERE
                istep = 11332
                iinfo = -1
                Goto 9999
             endif
             tempsepi(t+1)=jptr+1

             ! -------------------------------------
             ! [1.1.3.3.3] Generate weight of the separator
             !  to be sended put value on vtxpos
             ! -------------------------------------

             vtxpos(t) = vtxweight(vtxid)

          enddo

          nnzaintrf = jptr
          domintrfnnz(pcnb+1) = nnzaintrf

          !-------------------------------------------------------------
          ! [dbg.1]  Write datastructure to files
          !-------------------------------------------------------------

          debug=0
          if(debug.eq.1)then

             !------------------------------------------------
             ! [dbg.1.0] Initialize variables, etc.
             !------------------------------------------------

             nbsep = nproc - 1
             nbvi =m-1

             !------------------------------------------------
             ! [dbg.1.1]  Write nbvi, indexVI, ptrindexVI
             !------------------------------------------------

             write(unit=pass,FMT='(I5)') pcnb
             write(unit=filename,FMT='(A,"indexvi_",A)')&
                  trim(pbname),adjustl(pass)

             open(unit=18,file=filename,status="REPLACE",iostat=ios1)
             write(unit=18,FMT='(I16)') nbvi
             do p=1,nbvi
                write(unit=18,FMT='(I16)') indexVi(p)
             enddo
             do p=1,nbvi+1
                write(unit=18,FMT='(I16)') ptrindexVi(p)
             enddo
             close(unit=18)

             !------------------------------------------------
             ! [dbg.1.2]  Write IndexIntrf
             !------------------------------------------------

             write(unit=pass,FMT='(I5)') pcnb
             write(unit=filename,FMT='(A,"indexintrf_",A)')&
                  trim(pbname),adjustl(pass)
             open(unit=19,file=filename,status="REPLACE",iostat=ios1)
             do p=1,pt
                !       write(unit=19,FMT='(I16)') IndexIntrf(p)
                write(unit=19,FMT='(I16)') metperm(IndexIntrf(p))
             enddo
             close(unit=19)

             !------------------------------------------------
             ! [dbg.1.3]  Write My Interface
             !------------------------------------------------

             write(unit=pass,FMT='(I5)') pcnb
             write(unit=filename,FMT='(A,"myintrf_",A)')&
                  trim(pbname),adjustl(pass)
             open(unit=20,file=filename,status="REPLACE",iostat=ios1)
             do p=1,intrfsize
                !       write(unit=20,FMT='(I16)') listinterface(p)
                write(unit=20,FMT='(I16)') metperm(listinterface(p))
             enddo
             close(unit=20)

          endif

          !---------------------------------------------------
          ! [dbg.1.5]  Write data for MATLAB Testing
          !---------------------------------------------------

          debug=0
          if(debug.eq.1)then
             !   write interior indices to a file for matlab
             write(unit=pass,FMT='(I5)') pcnb+1
             write(unit=filename,FMT='(A,"_myintrf_",A)')&
                  trim(pbname),adjustl(pass)
             open(unit=20,file=filename,status="REPLACE",iostat=ios1)
             do p=1,intrfsize
                write(unit=20,FMT='(I16)') listinterface(p)
                ! write(unit=20,FMT='(I16)') metperm(listinterface(p))
             enddo
             close(unit=20)
          endif

          !-------------------------------------------------------------
          ! [1.2] Send the domain's structure
          !-------------------------------------------------------------

          nnzainterior = domintnnz(pcnb+1)
          nnzaintrf    = domintrfnnz(pcnb+1)
          !    print*,'SENDING TO DOMAIN NUMBER ',pcnb
          !&   ,nnzaintrf,nnzainterior,nnzaintrf+nnzainterior,
          !&   ' nnz from ',nnza

          ! --------------------------------------------------
          ! [1.2.1] Prepare general information about sending
          ! sending vectors and send vectors to proc pcnb
          ! -------------------------------------------------

          !> @par infoproc description :
          !!
          !! - (1) = size of interface of proc pcnb        
          !! - (2) = nb of neighbor of proc pcnb
          !! - (3) = length of Index_Intrfc of proc pcnb
          !!         should be equal to ptrindexVi(nbvi+1)-1
          !! - (4) = starting indicies of domain used to form local
          !!         permutation on each proc
          !! - (5) = UNUSED
          !! - (6) = length of the nnz of the separator row to be sended
          !!         to proc pcnb
          !! - (7) = length of the nnz of the interior row to be sended
          !!         to proc pcnb
          !! - (8) = length of the nnz of all my nnz (7+8)
          !! - (9) = length of total interface used to create gbtoloc
          !!         permutation
          !! - (10) = logical interface ndof of proc pcnb
          !! - (11) = interior ndof of proc pcnb
          !! - (12) = interface ndof of proc pcnb
          !! - (13) = allndof proc pcnb
          !! - (14) = allndof of the whole system (all proc)
          !! - (15) = number of the local interior Lagrange multiplier
          !! - (16) = symetry of the local matrix
          !!
          allsendnnz=allsendnnz+domintnnz(pcnb+1)+domintrfnnz(pcnb+1)
          infoproc(:) = -1

          nnzainterior = domintnnz(pcnb+1)
          nnzaintrf    = domintrfnnz(pcnb+1)

          infoproc(1) = intrfsize  
          infoproc(2) = m-1         
          infoproc(3) = pt          
          infoproc(4) = domstptr(pcnb+1)
          infoproc(6) = nnzaintrf          
          infoproc(7) = nnzainterior        
          infoproc(8) = nnzaintrf+nnzainterior
          infoproc(9) = totinterface        
          IF(RHSWAY.EQ.2) infoproc(10)= domlogicintrfdof(pcnb+1) 
          infoproc(11)= domintdof(pcnb+1)    
          infoproc(12)= domintrfdof(pcnb+1) 
          infoproc(13)= domintdof(pcnb+1) + domintrfdof(pcnb+1)
          infoproc(14)= ndof 
          infoproc(15)= domLg(pcnb+1)
          infoproc(16)= symtype


          if(pcnb.eq.me)then

             ! -----------------------------------------------
             ! [1.2.2] Master copy his domain's data
             ! -----------------------------------------------
             Do t=1,GINFOSIZE
                myinfoproc(t)  = infoproc(t)
             End Do
             mysizeIntrf    = myinfoproc(1)     ! size of interface of proc pcnb
             mynbvi         = myinfoproc(2)     ! nb of neighbor of proc pcnb
             lenindintrf    = myinfoproc(3)

             mynnzaintrf    = myinfoproc(6)     ! length of the nnz of the separator row to be sended to proc pcnb
             mynnzainterior = myinfoproc(7)     ! length of the nnz of the interior row to be sended to proc pcnb
             mynnza         = myinfoproc(8)     ! length of the nnz of all my nnz (7+8)
             gballintrf      = myinfoproc(9)    ! length of total interface used to create gbtoloc permutation
             myndoflogicintrf= myinfoproc(10)   ! logical interface ndof of proc pcnb
             myndofinterior = myinfoproc(11)    ! interior ndof of proc pcnb
             myndofintrf    = myinfoproc(12)    ! interface ndof of proc pcnb
             myndof         = myinfoproc(13)    ! allndof of proc pcnb
             gballndof      = myinfoproc(14)    ! allndof of the whole system (all proc)
             myint_Lg       = myinfoproc(15)    ! number of the local interior Lagrange multiplier
             symtype        = myinfoproc(16)

             ALLOCATE(myindexVi(mynbvi+1), stat=iinfo)
             if (iinfo /= 0 ) goto 959
             ALLOCATE(myptrindexVi(mynbvi+1), stat=iinfo)
             if (iinfo /= 0 ) goto 959
             ALLOCATE(myindexintrf(lenindintrf), stat=iinfo)
             if (iinfo /= 0 ) goto 959
             ALLOCATE(myinterface(mysizeIntrf+myint_Lg), stat=iinfo)
             if (iinfo /= 0 ) goto 959
             ALLOCATE(myiacsr(myndof+1), stat=iinfo)
             if (iinfo /= 0 ) goto 959
             BUGMUMP=mysizeIntrf+myint_Lg
             ALLOCATE(myja(mynnza+BUGMUMP), stat=iinfo)
             if (iinfo /= 0 ) goto 959
             ALLOCATE( mya(mynnza+BUGMUMP), stat=iinfo)
             if (iinfo /= 0 ) goto 959

             myinterface(:)=0
             myiacsr(:) = 0
             myja(:)    = 0
             mya(:)     = XMPH_FLOATZERO

             Call Icopy(mynbvi,indexVi,1,myindexVi,1)
             Call Icopy(mynbvi+1,ptrindexVi,1,myptrindexVi,1)
             Call Icopy(lenindintrf,IndexIntrf,1,myindexintrf,1)
             Call Icopy(mysizeIntrf,listinterface,1,myinterface,1)

             ! PAY ATTENTION that tempsepi use the same order 
             !               of indices as listinterface
             ed = myndofinterior
             myiacsr(ed+1)= mynnzainterior+1
             do t=2,mysizeIntrf+1
                rowsiz=tempsepi(t)-tempsepi(t-1)
                myiacsr(ed+t)= myiacsr(ed+t-1)+rowsiz
             enddo
             Call Icopy(mynnzaintrf,tempsepj,1,&
                  myja(mynnzainterior+1),1)
             Call XMPH_ARITHcopy(mynnzaintrf,tempsepval,1,&
                  mya(mynnzainterior+1),1)

          else

             ! -----------------------------------------------
             ! [1.2.3] Master send other domain's data
             ! -----------------------------------------------

             tagg = pcnb
             Call MPI_SEND(infoproc(1),GINFOSIZE,MPH_INTMPI,pcnb,&
                  tagg,comm,infompi)
             ASSRT(infompi == MPI_SUCCESS)

             Call MPI_SEND(indexVi(1),m-1,MPH_INTMPI,pcnb,&
                  tagg+1,comm,infompi)
             ASSRT(infompi == MPI_SUCCESS)

             Call MPI_SEND(ptrindexVi(1),m,MPH_INTMPI,pcnb,&
                  tagg+2,comm,infompi)
             ASSRT(infompi == MPI_SUCCESS)

             Call MPI_SEND(IndexIntrf(1),pt,MPH_INTMPI,pcnb,&
                  tagg+3,comm,infompi)
             ASSRT(infompi == MPI_SUCCESS)

             Call MPI_SEND(listinterface(1),intrfsize,MPH_INTMPI,pcnb,&
                  tagg+4,comm,infompi)
             ASSRT(infompi == MPI_SUCCESS)

             Call MPI_SEND(tempsepi(1),intrfsize+1,MPH_INTMPI,pcnb,&
                  tagg+10,comm,infompi)
             ASSRT(infompi == MPI_SUCCESS)

             Call MPI_SEND(tempsepj(1),nnzaintrf,MPH_INTMPI,pcnb,&
                  tagg+11,comm,infompi)
             ASSRT(infompi == MPI_SUCCESS)

             Call MPI_SEND(tempsepval(1),nnzaintrf,XMPH_FLOATMPI,&
                  pcnb,tagg+12,comm,infompi)
             ASSRT(infompi == MPI_SUCCESS)

          endif
       enddo

       ! -----------------------------------------------
       ! [1.2.4] Cleanup the sending step
       ! -----------------------------------------------

       DEALLOCATE(intrfindices)
       DEALLOCATE(intrfproc)
       DEALLOCATE(listinterface)
       DEALLOCATE(indexIntrf)
       DEALLOCATE(indexIntrfproc)
       DEALLOCATE(indexVi)
       DEALLOCATE(ptrindexVi)
       Deallocate(tempsepj)
       Deallocate(tempsepval)
       tdata=MPI_Wtime()-tdata

       IF((RHSWAY.EQ.1).AND.(scatpt.ne.combivtxpcsz)) goto 966 ! way1
       IF((RHSWAY.EQ.2).AND.(scatptlogic.ne.gballintrf)) goto 965 ! way2

    ELSE  ! I AM NOT THE MASTER

       !------------------------------------------------------------------------
       ! [1.3] Receive the domain's structure
       !------------------------------------------------------------------------

       !----------------------------------------------------------------
       ! [1.3.1] Receive the sizes of domain's structure
       !----------------------------------------------------------------
       tagg = me
       Call MPI_RECV(myinfoproc(1),GINFOSIZE, MPI_INTEGER,&
            0,tagg, comm,stat, infompi)

       mysizeIntrf   = myinfoproc(1)   ! size of interface
       mynbvi        = myinfoproc(2)   ! nb of neighbor
       lenindintrf   = myinfoproc(3)   ! length of Index_Intrfc
       domst         = myinfoproc(4)   ! starting indicies of my domain  in the global indicies.
       ! used to form local permutation on each proc

       mynnzaintrf    = myinfoproc(6)     ! length of the nnz of the separator row to be sended to proc pcnb
       mynnzainterior = myinfoproc(7)     ! length of the nnz of the interior row to be sended to proc pcnb
       mynnza         = myinfoproc(8)     ! length of the nnz of all my nnz (7+8)
       gballintrf      = myinfoproc(9)     ! length of total interface used to create gbtoloc permutation
       myndoflogicintrf= myinfoproc(10)    ! logical interface ndof of proc pcnb
       myndofinterior = myinfoproc(11)    ! interior ndof of proc pcnb
       myndofintrf    = myinfoproc(12)    ! interface ndof of proc pcnb
       myndof         = myinfoproc(13)    ! allndof of proc pcnb
       gballndof      = myinfoproc(14)    ! allndof of the whole system (all proc)
       myint_Lg       = myinfoproc(15)    ! number of the local interior Lagrange multiplier
       symtype        = myinfoproc(16)    ! symmetry of the local matrix

       !----------------------------------------------------------------
       ! [1.3.1] Allocate the domain's data
       !----------------------------------------------------------------

       ALLOCATE(myindexVi(mynbvi+1), stat=iinfo)
       if (iinfo /= 0 ) goto 959
       ALLOCATE(myptrindexVi(mynbvi+1), stat=iinfo)
       if (iinfo /= 0 ) goto 959
       ALLOCATE(myindexintrf(lenindintrf), stat=iinfo)
       if (iinfo /= 0 ) goto 959
       ALLOCATE(myinterface(mysizeIntrf+myint_Lg), stat=iinfo)
       if (iinfo /= 0 ) goto 959
       ALLOCATE(tempsepi(mysizeIntrf+1), stat=iinfo)
       if (iinfo /= 0 ) goto 959
       ALLOCATE(myiacsr(myndof+1), stat=iinfo)
       if (iinfo /= 0 ) goto 959
       BUGMUMP=mysizeIntrf+myint_Lg
       ALLOCATE(myja(mynnza+BUGMUMP), stat=iinfo)
       if (iinfo /= 0 ) goto 959
       ALLOCATE( mya(mynnza+BUGMUMP), stat=iinfo)
       if (iinfo /= 0 ) goto 959
       myinterface(:)=0
       myiacsr(:) = 0
       myja(:)    = 0
       mya(:)     = XMPH_FLOATZERO

       !----------------------------------------------------------------
       ! [1.3.2] Receive the domain's structure arrays
       !----------------------------------------------------------------

       CALL MPI_RECV(myindexVi(1),mynbvi, MPH_INTMPI,&
            0,tagg+1, comm,stat, infompi)
       ASSRT(infompi == MPI_SUCCESS)

       CALL MPI_RECV(myptrindexVi(1),mynbvi+1, MPH_INTMPI,&
            0,tagg+2, comm,stat, infompi)
       ASSRT(infompi == MPI_SUCCESS)

       CALL MPI_RECV(myindexintrf(1),lenindintrf,MPH_INTMPI,&
            0,tagg+3, comm,stat, infompi)
       ASSRT(infompi == MPI_SUCCESS)

       CALL MPI_RECV(myinterface(1),mysizeIntrf,MPH_INTMPI,&
            0,tagg+4, comm,stat, infompi)
       ASSRT(infompi == MPI_SUCCESS)

       CALL MPI_RECV(tempsepi(1),mysizeIntrf+1,MPH_INTMPI,&
            0,tagg+10,comm,stat, infompi)
       ASSRT(infompi == MPI_SUCCESS)

       CALL MPI_RECV(myja(mynnzainterior+1),mynnzaintrf,&
            MPH_INTMPI,0,tagg+11,comm,stat, infompi)
       ASSRT(infompi == MPI_SUCCESS)

       CALL MPI_RECV(mya(mynnzainterior+1),mynnzaintrf,&
            XMPH_FLOATMPI,0,tagg+12,comm,stat, infompi)
       ASSRT(infompi == MPI_SUCCESS)

       !----------------------------------------------------------------
       ! [1.3.3] Store tempsepi in iacsr(myndofinterior+1-->end)
       !----------------------------------------------------------------

       !----------------------------------------------------------------
       ! [1.3.3] Store tempsepi in iacsr(myndofinterior+1-->end)
       !----------------------------------------------------------------

       !   store tempsepi in iacsr(myndofinterior+1-->end)
       !   convert values in a way that the separator will be at the end
       !   that mean that separator indices start at mynnzainterior+1.
       !   NB: the iacsr(myndofinterior+1) value will be replaced
       !   during the scatterv by the received one but when
       !   converting the revceived value to the good ones
       !   starting from 1 its come back to mynnzainterior+1
       ! PAY ATTENTION that tempsepi use the same order of indices as listinterface

       !    ed = myndofinterior
       !    myiacsr(ed+1)= mynnzainterior+1
       !    do t=2,mysizeIntrf+1
       !       rowsiz=tempsepi(t)-tempsepi(t-1)
       !       myiacsr(ed+t)= myiacsr(ed+t-1)+rowsiz
       !    enddo

       ed = myndofinterior
       k = mynnzainterior ! amount to be added to continue numbering interface row after interior row
       do t=1,myndofintrf+1
          myiacsr(ed+t)= tempsepi(t)+ k
       enddo

       if(myiacsr(myndof+1).ne.(mynnza+1)) goto 969

    ENDIF  ! END IF I AM MASTER

    !---------------------------------------------------------------------------
    ! [1.4] Finish computing the domain structure, checks and prints
    !---------------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! [1.4.1] Test the sended nnz
    !-------------------------------------------------------------------

    istep = 141
    IF(ME.EQ.0) THEN
       ! TEST THE NNZ SENDED TO BE EQUAL TO NNZ INITIAL
       Call mph_logWithInfo(MSG_DEBUG, allsendnnz, &
       "Domain Decomposition -> total number of entries sended =")
       Call mph_logWithInfo(MSG_DEBUG, nnza, &
       "Domain Decomposition -> number of entries available =")

       ! if(rhsway.eq.2)print*,'domlogicintrfdof  :',(domlogicintrfdof(i),i=1,nproc)
       If(allsendnnz.ne.nnza)goto 960
    ENDIF

    !-------------------------------------------------------------------
    ! [1.4.2] Test lenindintrf
    !-------------------------------------------------------------------

    istep = 142
    if(lenindintrf.ne.myptrindexVi(mynbvi+1)-1) then
       iinfo = -1
       Goto 9999
    endif

    !-------------------------------------------------------------------
    ! [1.4.3] Inform user 
    !-------------------------------------------------------------------

    CALL MPI_BARRIER(comm,infompi)

    If (me.eq.0) Call mph_log(MSG_DEBUG, &
            "Domain Decomposition -> progress = data structure [DONE]")

    !---------------------------------------------------------------------------
    ! [2] Scatter the matrix interior entries
    !---------------------------------------------------------------------------

    tscat=MPI_Wtime()

    !----------------------------------------------------------------
    ! [2.1] Master Prepares the interior entries to be scattered
    !----------------------------------------------------------------

    ! Send matrix interior row [interface interior]
    IF(ME.EQ.0) THEN


       ALLOCATE(procintcnt (nproc), stat=iinfo)
       if (iinfo /= 0 ) goto 959
       ALLOCATE(procintdisp(nproc), stat=iinfo)
       if (iinfo /= 0 ) goto 959
       ALLOCATE(procintrfdisp(nproc), stat=iinfo)
       if (iinfo /= 0 ) goto 959
       ALLOCATE(procnzdisp (nproc), stat=iinfo)
       if (iinfo /= 0 ) goto 959

       !sendcount and displacement for iacsr
       procintcnt(:)  = 0
       procintdisp(:) = 0
       do domid=1,nproc
          procintcnt(domid) = domintdof(domid)+1 ! coz iacsr of ndof row is ndof+1
          if(domid.eq.1)then
             procintdisp(domid)=gballintrf ! Scatterv start at {displs[i]*extent(sendtype)} so no need to (+1)
          else
             procintdisp(domid)=procintdisp(domid-1)+&
                  procintcnt(domid-1)-1
          endif
       enddo

       !sendcount(allready done domintnnz) and displacement for (ja a)
       procnzdisp(:) = 0
       do domid=1,nproc
          if(domid.eq.1)then
             procnzdisp(domid)=nnzallinterface ! +1 ! Scatterv start at 0 so no need to (+1)
          else
             procnzdisp(domid)=procnzdisp(domid-1) + domintnnz(domid-1)
          endif
       enddo
    ELSE  ! I AM NOT THE MASTER

       !lty> obsolete ? bloc commented by Azzam
       !    ALLOCATE(domintnnz(1))
       !    ALLOCATE(procnzdisp(1))
       !    ALLOCATE( ja(1))
       !    domintnnz(:)  = 0
       !    procnzdisp(:) = 0
       !    ja(:)=0
       !lty<

       ! Allocate false aguments for MPI calls
       Allocate(iacsr(1)      )
       Allocate(ja   (1)      )
       Allocate( a   (1)      )
       Allocate(domintnnz(1)  )
       Allocate(procintcnt(1) )
       Allocate(procintdisp(1))
       Allocate(procnzdisp(1) )

       ! intialize
       iacsr       (:) = 0
       ja          (:) = 0
       a           (:) = 0
       domintnnz   (:) = 0
       procintcnt  (:) = 0
       procintdisp (:) = 0
       procnzdisp  (:) = 0

    ENDIF  ! END IF I AM MASTER


    !----------------------------------------------------------------
    ! [2.2] Scatter the rows    - iacsr (interior row)
    !----------------------------------------------------------------

    Call MPI_Scatterv(iacsr(1),procintcnt(1),procintdisp(1),&
         MPH_INTMPI,myiacsr(1),myndofinterior+1,&
         MPH_INTMPI,master,comm,ierr)

    !----------------------------------------------------------------
    ! [2.3] Scatter the columns - ja (interior row)
    !----------------------------------------------------------------

    Call MPI_Scatterv(ja(1),domintnnz(1),procnzdisp(1),&
         MPH_INTMPI,myja(1),mynnzainterior,&
         MPH_INTMPI,master,comm,ierr)

    !----------------------------------------------------------------
    ! [2.2] Scatter the values -  a (interior row)
    !----------------------------------------------------------------

    Call MPI_Scatterv(a(1),domintnnz(1),procnzdisp(1),&
         XMPH_FLOATMPI,mya(1),mynnzainterior,&
         XMPH_FLOATMPI,master,comm,ierr)

    !----------------------------------------------------------------
    ! [2.3] Inform user
    !----------------------------------------------------------------

    CALL MPI_BARRIER(comm,infompi)
    tscat  = MPI_Wtime()-tscat

    If (me.eq.0) Call mph_log(MSG_DEBUG, &
         "Domain Decomposition -> progress = scattering data [DONE]")

    !---------------------------------------------------------------------------
    ! [3] Convert to local indices 
    !---------------------------------------------------------------------------

    !
    !Received Matrix is ordered as follow
    ![ Abb  Abi]  separator_row(indices between 1<-->gballintrf)
    ![ Aib  Aii]  interior_row (indices from azz(i)%[rl]gst  --> azz(i)%[rl]ged
    !
    !I will convert order to
    !
    ![ Aii  Aib]  interior_row  indices from   1                   to    myndofinterior  (here are myndofinterior row)
    ![ Abi  Abb]  separator_row indices from   myndofinterior+1    to    myndof          (here are myndofintrf    row)
    !

    CALL MPI_BARRIER(comm,infompi)
    tloc   = MPI_Wtime()

    !----------------------------------------------------------------
    ! [3.1] Form myiacsr
    !----------------------------------------------------------------

    ! form myiacsr note that part of this myiacsr that
    ! correspond to interface was done during receiving data

    k=myiacsr(1)-1 ! amount to be subtracted to start at 1
    do i=1,myndofinterior+1
       myiacsr(i)= myiacsr(i)- k
    enddo
    if(myiacsr(myndofinterior+1).ne.(mynnzainterior+1)) goto 970

    !----------------------------------------------------------------
    ! [3.2] Form permutation
    !----------------------------------------------------------------

    ! first create gbtoloc for the interface indicies
    ! in a way that interface are at end.[interior interface]
    ! then it is easy for interior just a minus
    ! PAY ATTENTION that I should form gbtoloc IN THE SAME ORDER
    ! that tempsepi was formed

    ALLOCATE(gbtoloc(gballintrf), stat=iinfo)
    if (iinfo /= 0 ) goto 959

    gbtoloc(:)=-1
    do l=1,mysizeIntrf
       gbtoloc(myinterface(l)) = myndofinterior+l
    enddo

    !----------------------------------------------------------------
    ! [3.3] Convert indicies and weight value
    !----------------------------------------------------------------

    domst = myinfoproc(4)
    domed = myinfoproc(4)+myndofinterior-1
    k=domst-1
    do j=1,mynnza
       !   convert indices
       jvp=myja(j)
       if(jvp.gt.domed) goto 967 !test to be moved after coz impossible to have jvp bigger than domed

       if(jvp.ge.domst) then ! column corresponding to interior (column of Abi or Aii on received matrix)
          ! so schift them to left to the first
          myja(j)= jvp - k
       else                     ! column corresponding to interface (column of Aib or Abb on received matrix
          ! so schit them to right using gbtoloc.
          ! why? look upper? coz interface jcol are not successive
          if(jvp.gt.gballintrf) goto 967
          myja(j)= gbtoloc(jvp)
       endif
    enddo

    !convert data structure
    do w=1,mynbvi
       do u=myptrindexVi(w),myptrindexVi(w+1)-1
          myindexIntrf(u) = gbtoloc(myindexIntrf(u))- myndofinterior
       enddo
    enddo
    !convert weight
    ! weight is already done coz it is done as the
    ! same order as myinterface so by order of interface

    !---------------------------------------------------------------------------
    ! [4] Compute information that I need
    !---------------------------------------------------------------------------

    ALLOCATE(myia(mynnza+BUGMUMP), stat=iinfo)
    if (iinfo /= 0 ) goto 959

    do i=1,myndof
       do j=myiacsr(i),myiacsr(i+1)-1
          myia(j)=i
       enddo
    enddo
    !FOR MUMPS BUGGING
    k=myndof-BUGMUMP
    do i=mynnza+1,mynnza+BUGMUMP
       k=k+1
       myia(i)= k
       myja(i)= k
       mya(i) = XMPH_FLOATZERO
    enddo


    IF(RHSWAY.EQ.1)THEN
       !============================================
       !               |||| RHS |||| way1
       !============================================
       ! Master INFORMATION USED TO SEND RECEIVE RHS SOL
       !============================================
       !============================================
       !rhs in sended (sol is received) using 2 scatterv (2 gatherv)
       !interface RHS are groupped into a tempvector and scattered
       !interior RHS are scattered directly
       !I will use the sorted scatindices [intrftoproc0 intrftoproc1 intrftoproc2 ...]
       !to group interface RHS then domintrfdof procintrfdisp for sending using scatterv
       !scatindices   : will contain the indices of the interface of all proc and sorted by proc.
       !                [intrftoproc0 intrftoproc1 intrftoproc2 ...]
       !                PAY ATTENTION THAT FOR EACH PROC:
       !                - indices should be as the same order as numbered by the proc itself.
       !                - indices should stored after permutation. it is easy to scatterv
       !I will use domintdof procintdisp to scatterv interior
       !SAME for gathering the solution
       !============================================
       IF(ME.EQ.0) THEN
          !============================================
          !============================================
          !domintdof procintdisp for interior dof
          procintdisp(:) = 0
          do domid=1,nproc
             if(domid.eq.1)then
                procintdisp(domid)=gballintrf ! Scatterv start at {displs[i]*extent(sendtype)} so no need to (+1)
             else
                procintdisp(domid)=procintdisp(domid-1)+ domintdof(domid-1)
             endif
          enddo
          !domintrfdof procintrfdisp for interface dof
          procintrfdisp(:) = 0
          do domid=1,nproc
             if(domid.eq.1)then
                procintrfdisp(domid)= 0 ! Scatterv start at {displs[i]*extent(sendtype)} so no need to (+1)
             else
                procintrfdisp(domid)= procintrfdisp(domid-1)&
                     + domintrfdof(domid-1)
             endif
          enddo
          !============================================
       ENDIF  ! END IF I AM MASTER
       !============================================
       !               |||| END RHS ||||
       !============================================
    ENDIF !RHSWAY



    IF(RHSWAY.EQ.2)THEN
       !============================================
       !               |||| RHS |||| way2
       !============================================
       ! Master INFORMATION USED TO SEND RECEIVE RHS SOL
       !============================================
       !============================================
       !rhs in sended (sol is received) using 2 scatterv (2 gatherv)
       !LOGICAL interface RHS are groupped into a tempvector and scattered
       !interior RHS are scattered directly
       !I will use the sorted scatlogicindices [logicintrftoproc0 logicintrftoproc1 logicintrftoproc2 ...]
       !to group interface RHS then domlogicintrfdof procintrfdisp for sending using scatterv
       !scatlogicindices   : will contain the logical indices of the interface of all proc and sorted by proc.
       !                     [logicintrftoproc0 logicintrftoproc1 logicintrftoproc2 ...]
       !                     PAY ATTENTION THAT FOR EACH PROC:
       !                     - indices should be as the same order as numbered by the proc itself.
       !                     - indices should stored after permutation. it is easy to scatterv
       !I will use domintdof procintdisp to scatterv interior
       !SAME for gathering the solution BUT each proc should
       !generate itslogical interface indices (mylogicintrf)
       !to know which indices should gather.
       !============================================
       IF(ME.EQ.0) THEN
          !============================================
          !domintdof procintdisp for interior dof
          procintdisp(:) = 0
          do domid=1,nproc
             if(domid.eq.1)then
                procintdisp(domid)=gballintrf ! Scatterv start at {displs[i]*extent(sendtype)} so no need to (+1)
             else
                procintdisp(domid)=procintdisp(domid-1)+ domintdof(domid-1)
             endif
          enddo
          !domlogicintrfdof procintrfdisp for interface dof
          procintrfdisp(:) = 0
          do domid=1,nproc
             if(domid.eq.1)then
                procintrfdisp(domid)= 0 ! Scatterv start at {displs[i]*extent(sendtype)} so no need to (+1)
             else
                procintrfdisp(domid)= procintrfdisp(domid-1)&
                     + domlogicintrfdof(domid-1)
             endif
          enddo
          !*     testing proc 0 (domain 1) should not have any logical interface
          !*     it mean that domlogicintrfdof(1)=0 and procintrfdisp(1:2)=0
          ! if(procintrfdisp(2).ne.0)goto 963

          !testing proc 0 (domain 1) should have all its interface as logical
          !it mean that domlogicintrfdof(1)=sizeinterface.
          !also proc N (domain N) should not have any logical interface
          !it mean that domlogicintrfdof(nproc)=0 and procintrfdisp(nproc)=procintrfdisp(nproc-1)
          if(procintrfdisp(2).ne.domintrfdof(1))goto 961
          !============================================
       ENDIF  ! END IF I AM MASTER
       !============================================
       !
       !
       !! ---------------------------------------------
       !! generate mylogicintrf in SAME ORDER AS
       !! scatlogicindices for SENDING RECEIVING RHS SOL
       !! ---------------------------------------------
       !keep indices if I am the lowest procid that share it.
       !first put 2 on the interface indices that should belong to a processor lowest than me.
       !then find the remaining indices (that I should take them) IN SAME ORDER AS WAS DEFINED IN scatlogicindices .
       !PAY ATTENTION :
       ! - I should use the MY LOCAL NUMBERING BUT IN SAME ORDER
       !   AS WAS DEFINED IN scatlogicindices
       !   scatlogicindices is formed from (listinterface == myinterface).
       !   by the way, listinterface contains the indices in order so I can form
       !   mylogicintrf directly without need to templogic vector (tempsepi)
       !   but I keep it coz it is more general and to avoid any problem if
       !   someday I change the order of listinterface. by this way I keep
       !   all time same order without having information about listinterface.
       ! - myindexIntrf contains LOCAL numbering of interface indices
       !   but starting from 0 instead of ndofinterior+1.
       !   it is very good for economize memory => SO I should add (+ndofinterior)
       !   to the indices value when forming mylogicintrf.
       !
       !testing proc 0 (domain 1) should not have any logical interface
       if((me.eq.(nproc-1)).and.(myndoflogicintrf.ne.0)) goto 964
       !
       ALLOCATE(mylogicintrf(myndoflogicintrf),STAT=ierr1)
       if(ierr1.ne.0) goto 962
       Do w=1, myndoflogicintrf
          mylogicintrf(w)=0
       End Do
       !
       Do w=1, mysizeIntrf+1
          tempsepi(w)=0
       End Do
       !
       do w=1,mynbvi
          if(myindexVi(w).lt.me)then
             do u=myptrindexVi(w),myptrindexVi(w+1)-1
                vtxid = myindexIntrf(u)
                tempsepi(vtxid)=2
             enddo
          endif
       enddo
       !-------------------------------------
       !way2 to generate the gatherv indices
       !-------------------------------------
       scatptlogic = 0
       do w=1,mysizeIntrf
          vtxid = gbtoloc(myinterface(w))- myndofinterior ! gbtoloc give local numbering (schur at end)
          ! while indices in tempsepi are stored
          ! using (local numbering - myndofinterior)
          ! to economize memory so I have to substract myndofinterior
          if(tempsepi(vtxid).ne.2)then
             scatptlogic = scatptlogic +1
             mylogicintrf(scatptlogic) = gbtoloc(myinterface(w)) !local numbered (schur at the end)
          endif
       enddo

       !-------------------------------------
       !============================================
       !               |||| END RHS ||||
       !============================================
    ENDIF !RHSWAY
    !============================================
    !               |||| END RHS ||||
    !============================================

    tloc   = MPI_Wtime()-tloc
    tparti = MPI_Wtime()-tparti

    !---------------------------------------------------------------------------
    ! [5] Write Timings
    !---------------------------------------------------------------------------

    IF(ME.EQ.0) THEN
       ! Retrieve Previous Timings
       tmetis = ana_timing(1)
       ttmp   = ana_timing(2)
       tsort  = ana_timing(3)

       ! Save new Timings
       ana_timing(4) = tdata
       ana_timing(5) = tscat
       ana_timing(6) = tloc
       ana_timing(7) = tparti

       ! Inform user
       msg_class=MSG_VERBOSE2

       Write(msg,'(2A,F14.5)') "Domain Decomposition -> Timing -> ", &
            "total                        =", tparti
       Call mph_log(msg_class,msg)

       Write(msg,'(2A,F14.5)') "Domain Decomposition -> Timing -> ", &
            "ordering the matrix          =", tmetis
       Call mph_log(msg_class,msg)

       Write(msg,'(2A,F14.5)') "Domain Decomposition -> Timing -> ", &
            "computing DD from ordering   =", ttmp-tsort
       Call mph_log(msg_class,msg)

       Write(msg,'(2A,F14.5)') "Domain Decomposition -> Timing -> ", &
            "computing the scattered data =", tdata
       Call mph_log(msg_class,msg)

       Write(msg,'(2A,F14.5)') "Domain Decomposition -> Timing -> ", &
            "scattering the data          =", tscat
       Call mph_log(msg_class,msg)

       Write(msg,'(2A,F14.5)') "Domain Decomposition -> Timing -> ", &
            "computing the local data     =", tloc
       Call mph_log(msg_class,msg)

    ENDIF  ! END IF I AM MASTER


    !---------------------------------------------------------------------------
    ! [dbg.2] Write local matrice with different indices
    !---------------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! [dbg.2.1]  Write matrix to a file  global indices using iacsr
    !-------------------------------------------------------------------

    ![ Aii  Aib]
    ![ Abi  Abb]
    debug=0

    !---------------------------------------------------------
    !  [dbg.2.1.1] Real case
    !---------------------------------------------------------

#if  defined(MAPHYS_ARITH_d) || defined(MAPHYS_ARITH_s)

    if(debug.eq.1)then
       write(unit=pass,FMT='(I5)') me
       write(unit=filename,FMT='(A,"permutated_",A)')&
            trim(pbname),adjustl(pass)
       open(unit=unitnb5,file=filename,status="REPLACE",iostat=ios1)

       write(unit=unitnb5,FMT='(A)')&
            "%%MatrixMarket matrix coordinate REAL symmetric"
       write(unit=unitnb5,FMT='(I16,I16,I16)')&
            gballndof,gballndof,(myiacsr(myndof+1)-1)
       domst = myinfoproc(4)-1
       ! for interior_row [Aii Aib]
       do i=1,myndofinterior
          do j=myiacsr(i),(myiacsr(i+1)-1)
             if(myja(j).le.myndofinterior) then !corresponding to interior Aii
                write(unit=unitnb5,FMT=027)&
                     i+domst,myja(j)+domst,mya(j)
             else ! correspond to interior<-->interface Aib
                write(unit=unitnb5,FMT=027)&
                     i+domst,myinterface(myja(j)-myndofinterior),mya(j)
             endif
          enddo
       enddo
       ! for interface_row [Abi Abb]
       do i=myndofinterior+1,myndof
          do j=myiacsr(i),(myiacsr(i+1)-1)
             if(myja(j).le.myndofinterior) then !correspond to interior<-->interface Abi
                write(unit=unitnb5,FMT=027)&
                     myinterface(i-myndofinterior),myja(j)+domst,mya(j)
             else ! correspond to interface Abb
                write(unit=unitnb5,FMT=027)&
                     myinterface(i-myndofinterior),&
                     myinterface(myja(j)-myndofinterior),mya(j)
             endif
          enddo
       enddo
       close(unit=unitnb5)
    endif

#endif

    !---------------------------------------------------------
    !  [dbg.2.1.2] Complex case
    !---------------------------------------------------------

#if defined(MAPHYS_ARITH_z) || defined(MAPHYS_ARITH_c)

    if(debug.eq.1)then
       write(unit=pass,FMT='(I5)') me
       write(unit=filename,FMT='(A,"permutated_",A)')&
            trim(pbname),adjustl(pass)
       open(unit=unitnb5,file=filename,status="REPLACE",iostat=ios1)
       write(unit=unitnb5,FMT='(A)')&
            "%%MatrixMarket matrix coordinate COMPLEX symmetric"
       write(unit=unitnb5,FMT='(I16,I16,I16)')&
            gballndof,gballndof,(myiacsr(myndof+1)-1)
       domst = myinfoproc(4)-1
       ! for interior_row [Aii Aib]
       do i=1,myndofinterior
          do j=myiacsr(i),(myiacsr(i+1)-1)
             if(myja(j).le.myndofinterior) then !corresponding to interior Aii
                write(unit=unitnb5,FMT=028)&
                     i+domst,myja(j)+domst,REAL(mya(j)),AIMAG(mya(j))
             else ! correspond to interior<-->interface Aib
                write(unit=unitnb5,FMT=028)&
                     i+domst,myinterface(myja(j)-myndofinterior),&
                     REAL(mya(j)),AIMAG(mya(j))
             endif
          enddo
       enddo
       ! for interface_row [Abi Abb]
       do i=myndofinterior+1,myndof
          do j=myiacsr(i),(myiacsr(i+1)-1)
             if(myja(j).le.myndofinterior) then !correspond to interior<-->interface Abi
                write(unit=unitnb5,FMT=028)&
                     myinterface(i-myndofinterior),myja(j)+domst,&
                     REAL(mya(j)),AIMAG(mya(j))
             else ! correspond to interface Abb
                write(unit=unitnb5,FMT=028)&
                     myinterface(i-myndofinterior),&
                     myinterface(myja(j)-myndofinterior),&
                     REAL(mya(j)),AIMAG(mya(j))
             endif
          enddo
       enddo
       close(unit=unitnb5)
    endif

#endif

    !-------------------------------------------------------------------
    ! [dbg.2.2] Write matrix to a file global indices using myia
    !-------------------------------------------------------------------

    ![ Aii  Aib]
    ![ Abi  Abb]
    !
    debug=0
#if  defined(MAPHYS_ARITH_d) || defined(MAPHYS_ARITH_s)

    !---------------------------------------------------------
    !  [dbg.2.2.1] Real case
    !---------------------------------------------------------

    if(debug.eq.1)then
       write(unit=pass,FMT='(I5)') me
       write(unit=filename,FMT='(A,"permutated_",A)')&
            trim(pbname),adjustl(pass)
       open(unit=unitnb5,file=filename,status="REPLACE",iostat=ios1)
       write(unit=unitnb5,FMT='(A)')&
            "%%MatrixMarket matrix coordinate REAL general"
       write(unit=unitnb5,FMT='(I16,I16,I16)')&
            gballndof,gballndof,(myiacsr(myndof+1)-1)
       domst = myinfoproc(4)-1
       do i=1,mynnza
          if(myia(i).le.myndofinterior) then ! interior_row [Aii Aib]
             if(myja(i).le.myndofinterior) then !corresponding to interior Aii
                write(unit=unitnb5,FMT=027)&
                     myia(i)+domst,myja(i)+domst,mya(i)
             else ! correspond to interior<-->interface Aib
                write(unit=unitnb5,FMT=027)&
                     myia(i)+domst,myinterface(myja(i)-myndofinterior),&
                     mya(i)
             endif
          else ! interface_row [Abi Abb]
             if(myja(i).le.myndofinterior) then !correspond to interior<-->interface Abi
                write(unit=unitnb5,FMT=027)&
                     myinterface(myia(i)-myndofinterior),myja(i)+domst,&
                     mya(i)
             else ! correspond to interface Abb
                write(unit=unitnb5,FMT=027)&
                     myinterface(myia(i)-myndofinterior),&
                     myinterface(myja(i)-myndofinterior),&
                     mya(i)
             endif
          endif
       enddo
       close(unit=unitnb5)
    endif

#endif

    !---------------------------------------------------------
    !  [dbg.2.2.2] Complex case
    !---------------------------------------------------------

#if defined(MAPHYS_ARITH_z) || defined(MAPHYS_ARITH_c)

    if(debug.eq.1)then
       write(unit=pass,FMT='(I5)') me
       write(unit=filename,FMT='(A,"permutated_",A)')                   &
            trim(pbname),adjustl(pass)
       open(unit=unitnb5,file=filename,status="REPLACE",iostat=ios1)
       write(unit=unitnb5,FMT='(A)')                                    &
            "%%MatrixMarket matrix coordinate COMPLEX symmetric"
       write(unit=unitnb5,FMT='(I16,I16,I16)')                          &
            gballndof,gballndof,(myiacsr(myndof+1)-1)
       domst = myinfoproc(4)-1
       do i=1,mynnza
          if(myia(i).le.myndofinterior) then ! interior_row [Aii Aib]
             if(myja(i).le.myndofinterior) then !corresponding to interior Aii
                write(unit=unitnb5,FMT=028)                             &
                     myia(i)+domst,myja(i)+domst,                       &
                     REAL(mya(i)),AIMAG(mya(i))
             else ! correspond to interior<-->interface Aib
                write(unit=unitnb5,FMT=028)                             &
                     myia(i)+domst,myinterface(myja(i)-myndofinterior), &
                     REAL(mya(i)),AIMAG(mya(i))
             endif
          else ! interface_row [Abi Abb]
             if(myja(i).le.myndofinterior) then !correspond to interior<-->interface Abi
                write(unit=unitnb5,FMT=028)                             &
                     myinterface(myia(i)-myndofinterior),myja(i)+domst, &
                     REAL(mya(i)),AIMAG(mya(i))
             else ! correspond to interface Abb
                write(unit=unitnb5,FMT=028)                             &
                     myinterface(myia(i)-myndofinterior),               &
                     myinterface(myja(i)-myndofinterior),               &
                     REAL(mya(i)),AIMAG(mya(i))
             endif
          endif
       enddo
       close(unit=unitnb5)
    endif

#endif

    !---------------------------------------------------------------------------
    ! [dbg.3] Write matrix to a file local indices
    !---------------------------------------------------------------------------

    ![ Aii  Aib]
    ![ Abi  Abb]

    debug=0
    if(debug.eq.1)then
       write(unit=pass,FMT='(I5)') me
       write(unit=filename,FMT='(A,"permutated_",A)')&
            trim(pbname),adjustl(pass)
       open(unit=unitnb5,file=filename,status="REPLACE",iostat=ios1)

#if  defined(MAPHYS_ARITH_d) || defined(MAPHYS_ARITH_s)
       write(unit=unitnb5,FMT='(A)') &
            "%%MatrixMarket matrix coordinate REAL symmetric"
       write(unit=unitnb5,FMT='(I16,I16,I16,I16)')&
            myndof,myndof,mynnza,mysizeIntrf
       do i=1,mynnza
          write(unit=unitnb5,FMT=027) myia(i),myja(i),mya(i)
       enddo
#endif

#if defined(MAPHYS_ARITH_z) || defined(MAPHYS_ARITH_c)
       write(unit=unitnb5,FMT='(A)') &
            "%%MatrixMarket matrix coordinate COMPLEX symmetric"
       write(unit=unitnb5,FMT='(I16,I16,I16,I16)')&
            myndof,myndof,mynnza,mysizeIntrf
       do i=1,mynnza
          write(unit=unitnb5,FMT=028) myia(i),myja(i), &
               REAL(mya(i)),AIMAG(mya(i))
       enddo
#endif

       close(unit=unitnb5)
    endif


    !---------------------------------------------------------------------------
    ! [dbg.3] Write datastructure to files
    !---------------------------------------------------------------------------

    debug=0
    if(debug.eq.1)then

       !----------------------------------------------------------------
       ! [dbg.3.1] Write nbvi, indexVI, ptrindexVI
       !----------------------------------------------------------------

       write(unit=pass,FMT='(I5)') me
       write(unit=filename,FMT='(A,"indexvi_",A)')&
            trim(pbname),adjustl(pass)

       open(unit=unitnb5,file=filename,status="REPLACE",iostat=ios1)
       write(unit=unitnb5,FMT='(I16)') mynbvi
       do p=1,mynbvi
          write(unit=unitnb5,FMT='(I16)') myindexVi(p)
       enddo
       do p=1,mynbvi+1
          write(unit=unitnb5,FMT='(I16)') myptrindexVi(p)
       enddo
       close(unit=unitnb5)

       !----------------------------------------------------------------
       ! [dbg.3.2] Write IndexIntrf
       !----------------------------------------------------------------

       write(unit=pass,FMT='(I5)') me
       write(unit=filename,FMT='(A,"indexintrf_",A)')&
            trim(pbname),adjustl(pass)
       open(unit=unitnb5,file=filename,status="REPLACE",iostat=ios1)
       do p=1,lenindintrf
          write(unit=unitnb5,FMT='(I16)') myIndexIntrf(p)
       enddo
       close(unit=unitnb5)

       !----------------------------------------------------------------
       ! [dbg.3.3] Write my interface
       !----------------------------------------------------------------

       write(unit=pass,FMT='(I5)') me
       write(unit=filename,FMT='(A,"myintrf_",A)')&
            trim(pbname),adjustl(pass)
       open(unit=unitnb5,file=filename,status="REPLACE",iostat=ios1)
       do p=1,mysizeIntrf
          write(unit=unitnb5,FMT='(I16)') myinterface(p)
       enddo
       close(unit=unitnb5)

    endif

    CALL MPI_BARRIER(comm,infompi)
    If (me == 0) Call mph_log (MSG_DEBUG,&
         "Domain Decomposition -> progress = whole step [DONE]") 

    !---------------------------------------------------------------------------
    ! [6] Save computed data, free memory, etc. 
    !---------------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! [6.1] Free memory
    !-------------------------------------------------------------------
    ! for all processor
    Deallocate(myiacsr)
    Deallocate(tempsepi)
    If (Associated(listinterface)) Deallocate(listinterface)

    ! specific to a processor
    If (me.eq.0) then
       !    IF(metalgo.eq.3) DEALLOCATE(vwgt)
       !    DEALLOCATE(xadj)
       DEALLOCATE(adjncy)
       DEALLOCATE(domstptr)
       DEALLOCATE(vtxpos)
       DEALLOCATE(vtxweight)
       DEALLOCATE(domintnnz)
       DEALLOCATE(domintrfnnz)
       DEALLOCATE(belongingto)

    Else

       ! Deallocate false arguments in MPI calls

       Deallocate(iacsr         )
       Deallocate(ja            )
       Deallocate( a            )
       Deallocate(domintnnz     )
       Deallocate(procintcnt    )
       Deallocate(procintdisp   )

    Endif

    If (Associated(procintcnt)) Deallocate (procintcnt)
    If (Associated(procnzdisp)) Deallocate (procnzdisp)

    !-------------------------------------------------------------------
    ! [6.2] Save local matrix
    !-------------------------------------------------------------------


    Call XMPH_sm_ijv&
         (sm_local, myia,myja,mya,&
         myndof,myndof,mynnza, symtype, iinfo)

    !-------------------------------------------------------------------
    ! [6.3] Save the domain informations
    !-------------------------------------------------------------------

    one_domain%mynbvi       = mynbvi
    one_domain%lenindintrf  = lenindintrf
    one_domain%gballintrf   = gballintrf

    one_domain%myindexVi    => myindexVi
    one_domain%myptrindexVi => myptrindexVi
    one_domain%myindexintrf => myindexintrf
    one_domain%myinterface  => myinterface
    one_domain%mylogicintrf => mylogicintrf

    one_domain%myndof           = myndof
    one_domain%myndofinterior   = myndofinterior
    one_domain%myndofintrf      = myndofintrf
    one_domain%myndoflogicintrf = myndoflogicintrf
    one_domain%mysizeIntrf      = mysizeIntrf
    one_domain%myint_Lg         = myint_Lg

    !-------------------------------------------------------------------
    ! [6.4] Save RHS partitioning informations
    !-------------------------------------------------------------------

    part_rhs%gballndof        = gballndof
    part_rhs%gballintrf       = gballintrf
    part_rhs%gbtoloc      => gbtoloc
    part_rhs%nbdom = nproc
    if ( me == 0) then
       part_rhs%metperm  => metperm

       part_rhs%rhsway           =  rhsway
       part_rhs%combivtxpcsz     =  combivtxpcsz
       part_rhs%procintdisp      => procintdisp
       part_rhs%procintrfdisp    => procintrfdisp
       part_rhs%domLg            => domLg
       part_rhs%domintdof        => domintdof
       part_rhs%domintrfdof      => domintrfdof
       part_rhs%domlogicintrfdof => domlogicintrfdof
       part_rhs%scatindices      => scatindices
       part_rhs%scatlogicindices => scatlogicindices

    end if


    ! --------------------------------------------------------------------------
    ! [*.*] End
    ! --------------------------------------------------------------------------

9999 Continue

    ! Print error/warning messages

     If (iinfo /= 0) Then
        If ( iinfo > 0) msg_class = MSG_WARNING
        If ( iinfo < 0) msg_class = MSG_ERROR

        Select Case(istep) 
           !--------------------------------------
           ! errors/warnings based on istep.
           !--------------------------------------
        Case(100); write(msg,FMT='(I5,A,I5)')me,&
             ' : error in memory allocation, error = ',&
             iinfo
        Case(111); Write(msg,'(A,3I10)')&
             'ERROR the jth index is repeated for proc:', j,thisindex,thisproc
           Call mph_log(msg_class,msg)
           write(msg,'(A)')&
             '  ERROR 981: MATRIX ENTRIES ARE PROBABLY DUPLICATED COZ: ',&
             'an index belong to the same processor is founded twice or more'
        Case(11332); Write(msg,'(A)') &
             ' Azz_readpart ERROR : an interface index not belonging to anyone '
        Case(142); Write(msg,'(A)') &
             'Internal error on array prtrindexVi'

           !--------------------------------------
           ! Specific errors
           !--------------------------------------

        Case(959); write(msg,FMT='(I5,A,I5)')me,&
             ' : error in memory allocation, error = ',&
             iinfo

        Case(960); write(msg,FMT='(I5,A,2I10)')me,&
             '  Azz_const_struct ERROR 960 nnz sended and nnza are not equal',&
             allsendnnz,nnza

        Case(961); write(msg,FMT='(I5,A,I10)')me,&
             '  Azz_const_struct ERROR 961 PROC 0 HAS BAD LOGICAL INTERFACE',&
             myndoflogicintrf,domintrfdof(1)

        Case(962); write(msg,FMT='(I5,A,A,I10)')me,&
             '  Azz_const_struct ERROR 962 ERROR ON ALLOCATING mylogicintrf',&
             'look it is maybe proc 0 ',ierr1

        Case(964); write(msg,FMT='(I5,A,I10)')me,&
             '  Azz_const_struct ERROR 964 PROC N HAS LOGICAL INTERFACE',&
             myndoflogicintrf

        Case(965); write(msg,FMT='(I5,A,2I10)')me,&
             '  Azz_const_struct ERROR 965 scatptlogic # from gballintrf',&
             scatptlogic,gballintrf

        Case(966); write(msg,FMT='(I5,A,A,2I10)')me,&
             '  Azz_const_struct ERROR 966 scatpt # from combination of all',&
             ' interface ',scatpt,combivtxpcsz

        Case(967); write(msg,FMT='(I5,A,2I10)')me,&
             '  Azz_const_struct ERROR 967 ja bigger than domed ',&
             jvp,domed


        Case(969); write(msg,FMT='(I5,2A,2I10)')me,&
             '  Azz_const_struct ERROR 969 myiacsr(myndof+1) # ',&
             'from (mynnza+1)',&
             myiacsr(myndof+1),(mynnza+1)

        Case(970); write(msg,FMT='(I5,A,A,2I10)')me,&
             '  Azz_const_struct ERROR 970 myiacsr(myndofinterior+1) # ',&
             'from (mynnzainterior+1)',&
             myiacsr(myndofinterior+1),mynnzainterior+1

        Case(973); write(msg,FMT='(A60)')&
             'Azz_const_struct ERROR 973 a subdomain has no interface'

        Case(974); write(msg,FMT='(A60,I10,I10)')&
             'Azz_const_struct ERROR 974 st and ptrindexVi has ## values',&
             st,ptrindexVi(m-1)

        Case(975); write(msg,FMT='(A60,I10,I10)')&
             'Azz_const_struct ERROR 975 pt and ptrindexVi(nbvi+1)-1 has ##&
             &values',&
             pt,ptrindexVi(m)-1

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

959  istep=959 ; iinfo = -1; Goto 9999
960  istep=960 ; iinfo = -1; Goto 9999
961  istep=961 ; iinfo = -1; Goto 9999
962  istep=962 ; iinfo = -1; Goto 9999
964  istep=964 ; iinfo = -1; Goto 9999
965  istep=965 ; iinfo = -1; Goto 9999
966  istep=966 ; iinfo = -1; Goto 9999
967  istep=967 ; iinfo = -1; Goto 9999
969  istep=969 ; iinfo = -1; Goto 9999
970  istep=970 ; iinfo = -1; Goto 9999
973  istep=973 ; iinfo = -1; Goto 9999
974  istep=974 ; iinfo = -1; Goto 9999
975  istep=975 ; iinfo = -1; Goto 9999


     ! specific writing format.

#if  defined(MAPHYS_ARITH_d) || defined(MAPHYS_ARITH_s)
027 FORMAT(I10,I10,1PE26.18)
#endif
#if defined(MAPHYS_ARITH_z) || defined(MAPHYS_ARITH_c)
028 FORMAT(I10,I10,1PE26.18,1PE26.18)
#endif


  End Subroutine XMPH_PART_DistGlobalMatrix


  ! [+] routine : XMPH_Part_compute_interface_weight --------------------------------
  !
  !> Compute the weight of each node on the interface
  !! using the logic interface description.
  !!
  !!----
  !!
  !! @param [in,out] lc_domain 
  !!     The local domain description           
  !!     - On input , lc_domain contains valid : 
  !!         - myndoflogicintrf
  !!         - mylogicintrf
  !!         - mysizeintrf
  !!     - On output, append "lc_domain%weight"
  !!
  !!----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_Part_compute_interface_weight( lc_domain, info )

    !* Module(s) & co. *!
    Use MPH_part_type, Only : maphys_domain_t
    Use mph_log_mod
#if MAPHYS_DEBUG
      Use mph_dbg_mod
#endif
    Implicit None

    !* Arguments *!
    Type(maphys_domain_t), Intent(inout) :: lc_domain
    Integer              , Intent(  out) :: info 

    !* Local variables *!

    ! Scalars
    Integer                        :: iinfo  
    Integer                        :: istep
    Integer                        :: msg_class
    MPH_INT                     :: i
    MPH_INT                     :: sizeIntrf      
    MPH_INT                     :: ndofInterior
    MPH_INT                     :: ndofLogicIntrf 

    ! Arrays
    MPH_INT                   , Pointer :: LogicIntrf (:)     
    Real(kind=8)                 , Pointer :: weight     (:)
    
    ! Strings
    Character(len=MAPHYS_STRL)     :: rname = "XMPH_Part_compute_interface_weight"


    !- End of header------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Initialize variables
    !---------------------------------------------------------------------------
    
    !
    iinfo = 0
    Nullify(weight)    
    Nullify(LogicIntrf)

    ! Get local aliases
    LogicIntrf     => lc_domain%myLogicIntrf
    sizeIntrf      =  lc_domain%mySizeIntrf
    ndofLogicIntrf =  lc_domain%myNdofLogicIntrf
    ndofInterior   =  lc_domain%myNdofInterior

    ! Allocate weight
    istep = 11
    Allocate(weight(sizeIntrf), STAT=iinfo)
    If (iinfo > 0) iinfo = -iinfo
    If (iinfo < 0) Goto 9999
    
    !---------------------------------------------------------------------------
    ! [2] Compute weight
    !---------------------------------------------------------------------------

    Do i= 1,sizeIntrf
       weight(i) = 0.d0
    End Do

    Do i=1,ndoflogicintrf
       weight( LogicIntrf(i) - ndofInterior ) = 1.d0
    End Do

#if MAPHYS_DEBUG
    Call MPH_dbg_init()
    Write(dbg_unit,*) "rank           =", dbg_rank
    Write(dbg_unit,*) "sizeIntrf      =", sizeIntrf
    Write(dbg_unit,*) "ndoflogicintrf =", ndoflogicintrf
    Write(dbg_unit,*) "ndofInterior   =", ndofInterior
    Do i= 1,ndoflogicintrf
       Write(dbg_unit, *) logicintrf(i)
    End Do

    Do i= 1,sizeIntrf
       Write(dbg_unit, *) weight(i)
    End Do
#endif

    !---------------------------------------------------------------------------
    ! [3] Finish
    !---------------------------------------------------------------------------

    ! Save result
    lc_domain%weight => weight

    ! Print error/warnings

9999 Continue
    If ( iinfo /=  0 ) Then
       
       If ( iinfo > 0) msg_class = MSG_WARNING
       If ( iinfo < 0) msg_class = MSG_ERROR
       
       Select Case(istep) 
       Case(11); Call mph_logWithInfo (msg_class,sizeIntrf,Trim(rname)//&
            " failed to allocate weight, size =")
     End Select

  End If

  ! Set the return code

  If ( iinfo == 0 ) info =  0
  If ( iinfo <  0 ) info = -istep
  If ( iinfo >  0 ) info = +istep

  End Subroutine XMPH_Part_compute_interface_weight


End Module XMPH_Part_distmatrix_mod


