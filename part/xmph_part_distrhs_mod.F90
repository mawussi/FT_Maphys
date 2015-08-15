! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"
! [+] module : XMPH_part_distrhs_mod ------------------------------------------------
!
!> Module containing routines to distribute the second member.
Module XMPH_part_distrhs_mod

  !* List of routines
  Public :: XMPH_part_distrhs

  !* Private constants *!
  Character(len=MAPHYS_STRL), Private, Parameter :: FLNAME = &
       "xmph_part_distrhs_mod.F90"

Contains

  ! [+] routine : XMPH_PART_DistRhs --------------------------------------------
  !
  !> Distribute the centralized right-hand-side.
  !!
  !! Distribute the input centralized right-hand-side, into the domain 
  !! right-hand-side "lc_rhs".
  !!
  !! See doc_partitioning.inc for details.
  !!
  !!-----
  !!
  !! @param[in    ] part_rhs      data necessary to distribute gb_rhs
  !! @param[in    ] comm          the MPI communicator
  !! @param[in    ] lc_domain     local domain description
  !! @param[in    ] gb_rhs        the second member to distribute (size gballndof)
  !!
  !! @param[out   ] lc_rhs        the local rhs            (size myndof     ) 
  !! @param[out   ] info          the routine status
  !!
  !!-----
  !!
  !! @author Luc   Giraud
  !! @author Azzam Haidar
  !! @author Yohan Lee-tin-yien
  !!
  !! @version 0.1
  !!
  !! @todo Restructure Azzam's comments
  !! @todo Modify Outputs ?
  !! @todo Check allocation errors
  !!

  Subroutine XMPH_PART_DistRhs       ( & ! intents
       part_rhs, comm,            & ! in
       lc_domain, gb_rhs,         & ! in
       lc_rhs, info               & ! out
       )

    Use XMPH_dense_matrix_mod
    Use MPH_part_type
    Use mph_log_mod
    Use mph_error_mod
    Use XMPH_sparse_matrix_mod, Only : &
         XMPH_xset
    Implicit None
    Include 'mpif.h'

    !* Subroutine arguments *!
    Type(maphys_rhs_partition_t) , Intent(in   ) :: part_rhs
    Integer                      , Intent(in   ) :: comm
    Type(maphys_domain_t)        , Intent(in   ) :: lc_domain
    XMPH_FLOAT , Pointer       , Intent(in   ) :: gb_rhs (:)
    Type(XMPH_dense_matrix_t)         , Intent(  out) :: lc_rhs
    Integer                      , Intent(  out) :: info

    !* Local variables *!

    ! contants
    Integer, Parameter :: masterproc = 0

    ! scalars
    Integer    :: rank
    MPH_INT :: i
    MPH_INT :: gballintrf
    MPH_INT :: gballndof
    MPH_INT :: myndof,myndofinterior,myndoflogicintrf
    MPH_INT :: myndofintrf

    ! arrays
    MPH_INT   , pointer :: mylogicintrf(:)
    MPH_INT   , pointer :: metperm(:),scatlogicindices(:)
    MPH_INT   , pointer :: domintdof(:),domlogicintrfdof(:)
    MPH_INT   , pointer :: procintdisp(:),procintrfdisp(:)

    XMPH_FLOAT , pointer :: permutatedrhs (:)
    XMPH_FLOAT , pointer :: temprhsintrf  (:)
    XMPH_FLOAT , pointer :: templogicrhs  (:)

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Initialize local variables
    !---------------------------------------------------------------------------

    ! Pointers

    Nullify(templogicrhs )
    Nullify(permutatedrhs) 
    Nullify(temprhsintrf )
    Nullify(metperm)
    Nullify(scatlogicindices)     
    Nullify(domintdof)
    Nullify(domlogicintrfdof)
    Nullify(procintdisp)
    Nullify(procintrfdisp)

    ! Get MPI rank
    Call MPI_Comm_rank(comm,rank,info)
    ASSRTMPI(info)

    ! Link to part_rhs data
    If(rank == masterproc) Then
       gballintrf       =  part_rhs%gballintrf
       gballndof        =  part_rhs%gballndof

       metperm          => part_rhs%metperm         
       scatlogicindices => part_rhs%scatlogicindices     
       domintdof        => part_rhs%domintdof       
       domlogicintrfdof => part_rhs%domlogicintrfdof
       procintdisp      => part_rhs%procintdisp     
       procintrfdisp    => part_rhs%procintrfdisp   
    Else
       gballintrf       = 1
       gballndof        = 1
    End If

    ! Link to lc_domain data
    myndof           = lc_domain%myndof
    myndofinterior   = lc_domain%myndofinterior
    myndofintrf      = lc_domain%myndofintrf
    myndoflogicintrf = lc_domain%myndoflogicintrf

    mylogicintrf    => lc_domain%mylogicintrf

    ! Allocate / Initialize rhs

    !> @warning should "lc_rhs" be allocated before ?
    ! local rhs

    Call XMPH_DM_Nullify(lc_rhs, info)
    MPH_ONFAILURE_GOTO9999(info)

    Call XMPH_DM_Create(lc_rhs, myndof, 1, myndof, info)
    MPH_ONFAILURE_GOTO9999(info)

    Allocate(templogicrhs(Max(1, myndoflogicintrf)),stat= info)
    CHCKALLOC(info)
    MPH_ONFAILURE_GOTO9999(info)

    Call XMPH_xset(templogicrhs,1,myndoflogicintrf,XMPH_FLOATZERO)

    ! global rhs

    If( rank == masterproc ) Then

       Allocate(permutatedrhs(gballndof), stat = info) 
       CHCKALLOC(info)
       MPH_ONFAILURE_GOTO9999(info)
       Call XMPH_xset(permutatedrhs,1,gballndof,XMPH_FLOATZERO)

       Allocate(temprhsintrf(gballintrf), stat = info)
       CHCKALLOC(info)
       MPH_ONFAILURE_GOTO9999(info)
       Call XMPH_xset(temprhsintrf,1,gballintrf,XMPH_FLOATZERO)

    End If

    ! Allocate & Initialize dummy MPI arguments
    If( rank /= masterproc ) Then

       ! allocate 
       Allocate(&
            permutatedrhs   (1), &
            temprhsintrf    (1), &
            domlogicintrfdof(1), &
            domintdof       (1), &
            procintrfdisp   (1), &
            procintdisp     (1), STAT = info)
       CHCKALLOC(info)
       MPH_ONFAILURE_GOTO9999(info)

       ! give dummy values
       temprhsintrf    (1) = 0
       domlogicintrfdof(1) = 0
       domintdof       (1) = 0
       procintdisp     (1) = 0
       procintrfdisp   (1) = 0

    End If


    !--------------------------------------------------------------------------
    ! [2] Apply permutation
    !--------------------------------------------------------------------------
    ! Permute the rhs according to maphys permutation
    ! [allinterface interior0 interior1 ...]

    ! permutatedrhs=rhsglobal(metperm)

    If(rank == masterproc) Then

       Do i=1,gballndof
          permutatedrhs(i) = gb_rhs( metperm(i) )
       End Do

    End If

    !--------------------------------------------------------------------------
    ! [3] Extract interface rhs
    !--------------------------------------------------------------------------
    ! Extract the interface rhs and group them
    ! into temprhsintrf ordered by proc
    !
    If(rank == masterproc) Then

       Do i=1,gballintrf
          temprhsintrf(i) = permutatedrhs(scatlogicindices(i))
       End Do

    End If

    !---------------------------------------------------------------------------
    ! [4] Scatter the logical interface right-hand-side
    !---------------------------------------------------------------------------

    Call MPI_Scatterv(temprhsintrf,domlogicintrfdof,&
         procintrfdisp,XMPH_FLOATMPI,&
         templogicrhs,myndoflogicintrf,&
         XMPH_FLOATMPI,masterproc,comm,info)
    ASSRTMPI(info)

    Do i=1,myndoflogicintrf
       lc_rhs%v(mylogicintrf(i)) = templogicrhs(i)
    End Do

    !---------------------------------------------------------------------------
    ! [5] Scatter the interior right-hand-side
    !---------------------------------------------------------------------------

    Call MPI_Scatterv(permutatedrhs,domintdof,&
         procintdisp,XMPH_FLOATMPI,&
         lc_rhs%v,myndofinterior,&
         XMPH_FLOATMPI,masterproc,comm,info)
    ASSRTMPI(info)

    !---------------------------------------------------------------------------
    ! [6] Exit Routine
    !---------------------------------------------------------------------------

9999 Continue
    
    ! Free memory
    If(Associated(permutatedrhs   )) DeAllocate(permutatedrhs   )
    If(Associated(temprhsintrf    )) DeAllocate(temprhsintrf    )
    If(Associated(templogicrhs    )) DeAllocate(templogicrhs    )
    If( rank /= masterproc ) Then
       If(Associated(domlogicintrfdof)) DeAllocate(domlogicintrfdof)
       If(Associated(domintdof       )) DeAllocate(domintdof       )
       If(Associated(procintrfdisp   )) DeAllocate(procintrfdisp   )
       If(Associated(procintdisp     )) DeAllocate(procintdisp     )
    End If


  End Subroutine XMPH_PART_DistRhs

End Module XMPH_part_distrhs_mod



