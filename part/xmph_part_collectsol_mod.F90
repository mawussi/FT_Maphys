! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"

! [+] module : XMPH_part_collectsol_mod ------------------------------------------------
!
!> Module containing the routines to collect the solution.
Module XMPH_part_collectsol_mod

Contains


  ! [+] routine : XMPH_PART_CollectSol -------------------------------------------------
  !
  !> Collect the solution which is centralized on the host.
  !! 
  !!
  !! @param[in    ]  comm        the MPI Communicator used for communications
  !! @param[in,out]  part_rhs    informations about rhs partitionning
  !! @param[in    ]  lc_domain   the domain description 
  !! @param[in    ]  domain_sol  local solution to be gathered 
  !! @param[   out]  gb_sol      global solution
  !! @param[   out]  info        routine status
  !!
  !! @author Luc   Giraud
  !! @author Azzam Haidar
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_PART_CollectSol ( & ! intents
       part_rhs,               & ! inout
       comm,  lc_domain,       & ! in
       domain_sol,             &              
       gb_sol, info            & ! out
       )
    !* Modules *!

    Use MPH_part_type
    Use XMPH_dense_matrix_mod, Only : &
         XMPH_dense_matrix_t, & ! type
         XMPH_dm_create,    & ! routines
         XMPH_dm_nullify
    Use mph_log_mod
    Use mph_error_mod
#if MAPHYS_DEBUG
    Use mph_dbg_mod
#endif

    Implicit None
    Include "mpif.h"

    !* Arguments *!

    Type(maphys_rhs_partition_t) , Intent(inout) :: part_rhs
    Integer                      , Intent(in   ) :: comm
    Type(maphys_domain_t)        , Intent(in   ) :: lc_domain
    Type(XMPH_dense_matrix_t)         , Intent(in   ) :: domain_sol 
    Type(XMPH_dense_matrix_t)         , Intent(out  ) :: gb_sol 
    Integer                      , Intent(out  ) :: info

    !* Local variables *!

    ! scalars
    Integer :: master
    Integer ::  istep
    Integer ::  iinfo
    Integer ::  msg_class
    Integer ::  rank

    MPH_INT :: i
    ! new
    MPH_INT :: intrf_ndof
    MPH_INT :: interior_ndof
    MPH_INT :: gb_ndof
    MPH_INT :: lc_intrf_contrib_ndof

    ! Arrays
    MPH_INT   , Pointer ::   domintdof        (:)
    MPH_INT   , Pointer ::   domlogicintrfdof (:)
    MPH_INT   , Pointer ::   procintdisp      (:)
    MPH_INT   , Pointer ::   procintrfdisp    (:)

    ! new
    MPH_INT   , Pointer ::   gb_matrix_perm    (:)
    MPH_INT   , Pointer ::   lc_intrf_filter   (:)
    MPH_INT   , Pointer ::   intrf2gb_numbering(:)  

    XMPH_FLOAT , Pointer ::   intrf_sol         (:)
    XMPH_FLOAT , Pointer ::   gb_sol_permutated (:)
    XMPH_FLOAT , Pointer ::   lc_intrf_contrib  (:)


    ! Strings
    Character(len=MAPHYS_STRL)   :: rname = "PART_Collect_sol"
    !- End of header -------------------------------------------------------------

    !-----------------------------------------------------------------------------
    ! [1] Initialize local variables
    !-----------------------------------------------------------------------------
#if MAPHYS_DEBUG
    Call MPH_dbg_init()
    Call MPH_dbg_prefix_prepend(".")
    Call MPH_dbg_prefix_prepend(rname)
#endif   

    ! Nullify pointers

    Nullify(domintdof         )
    Nullify(domlogicintrfdof  )
    Nullify(procintdisp       )
    Nullify(procintrfdisp     )

    Nullify(gb_matrix_perm    )
    Nullify(lc_intrf_filter   )
    Nullify(intrf2gb_numbering)  

    Nullify(lc_intrf_contrib )
    Nullify(intrf_sol        )
    Nullify(gb_sol_permutated)
    
    !-------------------------------------------------------------------
    ! [1.1] Local Scalars & aliases
    !-------------------------------------------------------------------

    ! scalars
    master = 0

    ! MPI
    istep = 11
    CALL MPI_COMM_RANK(comm,rank,iinfo)
    Call mph_error_checkMPI(iinfo)
    If (iinfo < 0 ) Goto 9999

    ! sizes
    interior_ndof         = lc_domain%myndofinterior
    lc_intrf_contrib_ndof = lc_domain%myndoflogicintrf

    ! aliases (original)
    If ( rank == master ) Then
       intrf_ndof            = part_rhs%gballintrf
       gb_ndof               = part_rhs%gballndof
       procintdisp      => part_rhs%procintdisp     
       procintrfdisp    => part_rhs%procintrfdisp   
       domintdof        => part_rhs%domintdof       
       domlogicintrfdof => part_rhs%domlogicintrfdof
    Else
       intrf_ndof            = 1
       gb_ndof               = 1
       Allocate( domintdof       (1))        
       Allocate( procintrfdisp   (1))
       Allocate( procintdisp     (1))
       Allocate( domlogicintrfdof(1))
    End If

    ! aliases (renamed)
    lc_intrf_filter    => lc_domain%mylogicintrf
    If ( rank == master ) Then
       gb_matrix_perm     => part_rhs%metperm         
       intrf2gb_numbering => part_rhs%scatlogicindices
    End If

    !-------------------------------------------------------------------
    ! [1.2] Allocate memory
    !-------------------------------------------------------------------

    ! local contribution to the interface
    istep = 121
    Allocate ( lc_intrf_contrib(Max(1,lc_intrf_contrib_ndof)), stat= iinfo )
    If (iinfo > 0 ) iinfo = -iinfo
    If( iinfo <  0) Goto 9999

    ! global interface solution
    istep = 122
    Allocate ( intrf_sol (intrf_ndof), stat= iinfo )
    If (iinfo > 0 ) iinfo = -iinfo
    If( iinfo <  0) Goto 9999

    ! global permutated solution
    istep = 123
    Allocate ( gb_sol_permutated (gb_ndof), stat= iinfo )
    If (iinfo > 0 ) iinfo = -iinfo
    If( iinfo <  0) Goto 9999

    ! default values
    Do i=1, lc_intrf_contrib_ndof
       lc_intrf_contrib  (i) = XMPH_FLOATZERO
    End Do
    Do i=1, intrf_ndof
       intrf_sol         (i) = XMPH_FLOATZERO
    End Do
    Do i=1, gb_ndof
       gb_sol_permutated (i) = XMPH_FLOATZERO
    End Do

    ! global solution
    istep = 124
    If ( rank == master ) Then
       Call XMPH_DM_Create( gb_sol, gb_ndof, 1, gb_ndof, iinfo)
       If( iinfo < 0) Goto 9999
    Else
       Call XMPH_DM_Nullify( gb_sol, iinfo )
       If( iinfo < 0) Goto 9999
    End If

    !-----------------------------------------------------------------------------
    ! [2] Get the solution on the interface
    !-----------------------------------------------------------------------------

    ! Filter the domain solution to obtain the local contribution to the interface
#if MAPHYS_DEBUG
    Call MPH_dbg_set_file("lc_intrf_filter")
#endif   

    Do i = 1, lc_intrf_contrib_ndof
       lc_intrf_contrib( i ) = domain_sol%v( lc_intrf_filter(i) )
#if MAPHYS_DEBUG
       write(dbg_unit,*) i, 1, lc_intrf_filter(i)
#endif   

    End Do

    ! Gather all local contribution to form the global interface solution
    istep = 22

    Call MPI_Gatherv(&
         lc_intrf_contrib(1),lc_intrf_contrib_ndof, XMPH_FLOATMPI,&
         intrf_sol(1) ,domlogicintrfdof(1),procintrfdisp(1) , XMPH_FLOATMPI,&
         master,comm,iinfo) 
    Call mph_error_checkMPI(iinfo)
    If( iinfo < 0) Goto 9999

#if MAPHYS_DEBUG
    Call MPH_dbg_set_file("intrf_sol")
    Do i=1,intrf_ndof
       write(dbg_unit,*) i, 1,  intrf_sol(i)
    End Do
#endif   


    !-----------------------------------------------------------------------------
    ! [3] Get the solution on the interiors
    !-----------------------------------------------------------------------------

    istep = 23

    Call MPI_Gatherv(&
         domain_sol%v(1)     ,interior_ndof              ,XMPH_FLOATMPI,&
         gb_sol_permutated(1),domintdof(1),procintdisp(1),XMPH_FLOATMPI,&
         master,comm,iinfo)
    Call mph_error_checkMPI(iinfo)
    If( iinfo < 0) Goto 9999


#if MAPHYS_DEBUG
    If ( rank == master ) Then 
       Call MPH_dbg_set_file("interior_sol")
       Do i=1,gb_ndof
          write(dbg_unit,*) i, 1,  gb_sol_permutated(i)
       End Do
    End If
#endif   

    !-----------------------------------------------------------------------------
    ! [4] Master Forms the global permutated solution "sol_permutated" 
    !-----------------------------------------------------------------------------

    master_form_sol : If ( rank == master ) Then 

       ! Empty the solution on the interface 
       Do i=1,intrf_ndof
          gb_sol_permutated(i) = XMPH_FLOATZERO
       End Do

       ! Give the real values
       Do i=1,intrf_ndof
          gb_sol_permutated( intrf2gb_numbering(i) ) = intrf_sol(i)
       End Do

#if MAPHYS_DEBUG
       Call MPH_dbg_set_file("gb_sol_permutated")
       Do i=1,gb_ndof
          write(dbg_unit,*) i, 1,  gb_sol_permutated(i)
       End Do
#endif   

    End If master_form_sol


    !-----------------------------------------------------------------------------
    ! [5] Master permutes back "sol" to form the global solution "gb_sol"
    !-----------------------------------------------------------------------------
    !azz>
    ! Now permutatedrhs is according to my permutation that is:
    ! [allinterface interior0 interior1 ...]
    ! So I can get back the original solution using metperm
    !azz<

    master_permute_back : If ( rank == master ) Then 
       do i=1, gb_sol%m
          gb_sol%v(gb_matrix_perm (i)) = gb_sol_permutated(i)
       enddo
    End If master_permute_back

    !-----------------------------------------------------------------
    ! [6] Exit routine
    !-----------------------------------------------------------------

#if MAPHYS_DEBUG
    Call MPH_dbg_exit()
#endif

    ! Print error/warning messages
9999 Continue
    If ( iinfo /=  0 ) Then
       
       If ( iinfo > 0) msg_class = MSG_WARNING
       If ( iinfo < 0) msg_class = MSG_ERROR
       
       Select Case(istep) 
       Case(121,122,123)
          Call mph_logWithInfo (msg_class,istep,Trim(rname)//&
               " in an allocation at internal step = ")
       Case(11,22,23)
          Call mph_logWithInfo (msg_class,istep,Trim(rname)//&
               " in an MPI Call at internal step = ")
          Call mph_abort
       Case(124)
          Call mph_logWithInfo (msg_class,istep,Trim(rname)//&
               " in an call to dense matrix module at internal step = ")
       Case Default
          Call mph_logWithInfo (msg_class,istep, Trim(rname)//&
               " at internal step =")
       End Select
    End If

    ! set return code
    If ( iinfo == 0 ) info =  0
    If ( iinfo <  0 ) info = -istep
    If ( iinfo >  0 ) info = +istep

    ! Free memory 
    If ( Associated(lc_intrf_contrib )) Deallocate ( lc_intrf_contrib )
    If ( Associated(intrf_sol        )) Deallocate ( intrf_sol )
    If ( Associated(gb_sol_permutated)) Deallocate ( gb_sol_permutated )
    If ( rank /= 0)Then
       If ( Associated( domintdof       )) Deallocate( domintdof       )        
       If ( Associated( domlogicintrfdof)) Deallocate( domlogicintrfdof)
       If ( Associated( procintrfdisp   )) Deallocate( procintrfdisp   )
       If ( Associated( procintdisp     )) Deallocate( procintdisp     )
    End If

  End Subroutine XMPH_PART_CollectSol
End Module XMPH_part_collectsol_mod

