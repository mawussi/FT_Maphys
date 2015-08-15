! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"

  !> handle the distributed system
  Module XMPH_distsystem_mod

    !* Module(s) *!
    Use XMPH_maphys_type
    Use mph_error_mod

    !* No implicit typing *!
    Implicit None

    !* Defined routine(s) *!
    Public :: XMPH_distsystem_read
    Public :: XMPH_distsystem_buildSubBlocs
    Public :: XMPH_distsystem_compute

    !* Private constants *!
    Character(len=MAPHYS_STRL), Private, Parameter :: &
         FLNAME= "XMPH_ARITHmph_maphys_mod.F90"

  Contains

    !> read the distributed system from files
    !! 
    !! @param[in,out] mphs the maphys instance
    !!          on output, lc_domain,rhs_part, and sls_domain%sm_A
    !!                     are set + their related statistics
    !! @param[in]     domunit the unit to local domain file
    !! @param[in]     smunit  the unit to local matrix file
    !! @param[in]     rpunit  the unit to rhs partitioning file
    Subroutine XMPH_distsystem_read(&
         mphs, domunit, smunit, rpunit, info )
      
      !* Modules *!
      Use MPH_maphys_enum
      Use MPH_domain_mod
      Use MPH_rhs_partition_mod
      Use MPH_mem_mod
      Use XMPH_maphys_type
      Use XMPH_sparse_matrix_mod

      !* Arguments *!
      Type(XMPH_maphys_t), Intent(inout) :: mphs
      Integer, Intent(in) :: domunit
      Integer, Intent(in) :: smunit
      Integer, Intent(in) :: rpunit

      Integer, Intent(out) :: info

      !- End of header ------------------------------------------------------

      ! [-] Init

      mphs%rinfo (RINFO_TIMING_In2LcSys) =  0.d0
      mphs%iinfo (IINFO_STRAT_PART     ) = -1

      Call MPH_domain_nullify(mphs%lc_domain)
      Call MPH_rhspart_null(mphs%part_rhs)
      Call XMPH_sm_nullify(mphs%sls_domain%sm_A,info)
      MPH_ONFAILURE_RETURN(info)

      ! [-] Read

      Call MPH_domain_read(mphs%lc_domain,domunit,info)
      MPH_ONFAILURE_RETURN(info)

      Call XMPH_sm_mmread(mphs%sls_domain%sm_A,smunit,info)
      MPH_ONFAILURE_RETURN(info)

      If (rpunit < 0) Then 
         info = 1
         Return
      End If

      Call MPH_rhspart_read(mphs%part_rhs,rpunit,info)
      MPH_ONFAILURE_RETURN(info)

      ! [-] Check

      CHCKASSRT(mphs%lc_domain%gballintrf == mphs%part_rhs%gballintrf, info )
      MPH_ONFAILURE_RETURN(info)

      ! [-] Stats
      Call MPH_mem_add2usage(mphs%mem,&
           XMPH_sm_sizeof(mphs%sls_domain%sm_A))

    End Subroutine XMPH_distsystem_read
    !
    !---------------------------------------------------------------------------
    !
    !> build a the distributed system from 
    !! what was provided by the user.
    !!
    !! @param [in,out] mphs
    !!     - on input, mphs%domain is partially set
    !!                 mphs%sls_domain%sm_A is set
    !!     - on output, mphs%domain is fully set
    !!                  mphs%sm_A{i,b}{i,b} are set
    Subroutine XMPH_distsystem_buildSubBlocs(mphs, info)

      !* Modules *!
      Use MPH_mem_mod
      Use MPH_maphys_enum
      Use XMPH_maphys_type
      Use XMPH_sparse_matrix_mod

      !* Arguments *!
      Type(XMPH_maphys_t) , Intent(inout) :: mphs
      Integer, Intent(out) :: info

      !* External routine *!
      Real(Kind=8), External :: MPI_Wtime

      !* Local variables *!
      Real(Kind=8) :: timing
      MPH_INT :: Aii_ndof

      !- End of header ------------------------------------------------------

      ! [-] Init

      info = 0
      Aii_ndof = mphs%lc_domain%myndofinterior

      ! [-] Build the subblocs of the local matrix

      timing = MPI_Wtime()

      Call XMPH_sm_get_submatrices( & 
           Aii_ndof, mphs%sls_domain%sm_A, & 
           mphs%sm_Aii, mphs%sm_Aib,       & 
           mphs%sm_Abi, mphs%sm_Abb,       &
           info )
      timing = MPI_Wtime() - timing

      MPH_ONFAILURE_RETURN(info)

      ! [-] update counters

      mphs%rinfo (RINFO_TIMING_ExtBlocs) = timing

      Call MPH_mem_add2usage(mphs%mem,XMPH_sm_sizeof(mphs%sm_Aii))
      Call MPH_mem_add2usage(mphs%mem,XMPH_sm_sizeof(mphs%sm_Aib))
      Call MPH_mem_add2usage(mphs%mem,XMPH_sm_sizeof(mphs%sm_Abi))
      Call MPH_mem_add2usage(mphs%mem,XMPH_sm_sizeof(mphs%sm_Abb))

      mphs%iinfo (IINFO_LCMAT_SIZEOF  ) = Byte2MByte( &
           XMPH_sm_sizeof( mphs%sls_domain%sm_A ) + &
           XMPH_sm_sizeof( mphs%sm_Aii ) + &
           XMPH_sm_sizeof( mphs%sm_Aib ) + &
           XMPH_sm_sizeof( mphs%sm_Abi ) + &
           XMPH_sm_sizeof( mphs%sm_Abb ) )

      mphs%iinfo (IINFO_SCHUR_ORDER    ) = mphs%sm_Abb%m
      mphs%iinfo (IINFO_LCMAT_ORDER    ) = mphs%sls_domain%sm_A%m
      mphs%iinfo (IINFO_LCMAT_NBENTRIES) = mphs%sls_domain%sm_A%nnz

    End Subroutine XMPH_distsystem_buildSubBlocs
    !
    !---------------------------------------------------------------------------
    !
    !> Compute the distributed system from a global linear system (matrix)
    !!
    !!  distribute the global linear system into local linear systems (one
    !!  per MPI process).
    !!
    !!----
    !!
    !!  @param[in,out] mphs the maphys instance
    !!  @param[in,out] gb_A the global matrix to distribute   
    !!
    !!----
    !!
    !! @author Yohan Lee-tin-yien
    !!
    Subroutine XMPH_distsystem_compute(mphs,master,gb_A,comm,info)

      !* Modules & co. *!

      Use MPH_maphys_enum
      Use XMPH_maphys_type
      Use XMPH_sparse_matrix_mod
      Use MPH_part_type
      Use XMPH_part_mod
      Use MPH_mem_mod
      Implicit None
      Include 'mpif.h'

      !* Subroutine arguments *!
      Type(XMPH_maphys_t)       , Intent(inout) :: mphs
      Integer                   , Intent(in)    :: master
      Type(XMPH_sparse_matrix_t), Intent(inout) :: gb_A
      Integer                   , Intent(inout) :: comm
      Integer                   , Intent(  out) :: info

      !* Local Variables *!

      ! Scalars
      Integer :: partstrat
      Integer :: rank
      Integer :: nbdom ! the number of domains (= subsystems)
      Real(kind=8) :: time_Input2LocalSystem

      ! Derived types

      real(kind=8) :: ana_timing(ANA_TIMING_SIZE)

      Type(maphys_matrix_graph_t) :: graph_global_A
      Type(maphys_binary_tree_t)  :: part_tree
      Type(maphys_domains_t)      :: domains

      !- End of header----------------------------------------------------------

      !-------------------------------------------------------------------------
      ! [-] Init
      !-------------------------------------------------------------------------

      ! Scalars
      info = 0
      partstrat = mphs%ikeep (IKEEP_PART_Strategy)

      Call MPI_Comm_rank(mphs%comm,rank,info)
      ASSRT( info == MPI_SUCCESS )

      Call MPI_Comm_size(mphs%comm,nbdom,info)
      ASSRT( info == MPI_SUCCESS )

      ! Derived types
      Allocate(mphs%sls_domain%sm_A,STAT=info)
      Call XMPH_sm_nullify(mphs%sls_domain%sm_A, info)
      MPH_ONFAILURE_RETURN(info)

      !-------------------------------------------------------------------------
      ! [-] Distribute the system
      !-------------------------------------------------------------------------

      time_input2localsystem = MPI_Wtime()
      If (rank == master ) Then

         !----------------------------------------------------------------------
         ! [--] Order the system ( Modified METIS or SCOTCH )
         !----------------------------------------------------------------------

         Call XMPH_PART_OrderGlobalMatrix(          & ! intents
              nbdom, partstrat,           & ! in
              gb_A, ana_timing,                & ! inout 
              graph_global_A, part_tree, info & ! out    
              )
         MPH_ONFAILURE_ABORT(info)

         !----------------------------------------------------------------------
         ! [--] Change data representation, construct the boundaries of domains
         !----------------------------------------------------------------------

         Call XMPH_PART_Build_Domains(      & ! intents
              gb_A, graph_global_A,          & ! inout
              part_tree, ana_timing,         &
              domains, info                 & ! out
              ) 
         MPH_ONFAILURE_ABORT(info)

      endif ! rank == master

      !-------------------------------------------------------------------------
      ! [--] Distribute the domains, save necessary data to distribute the RHS
      !-------------------------------------------------------------------------

      Call XMPH_PART_DistGlobalMatrix(     & ! intents
           rank, comm, domains,              & ! in
           gb_A, graph_global_A, ana_timing, & ! inout
           mphs%part_rhs,                    & ! out 
           mphs%lc_domain, mphs%sls_domain%sm_A, &
           info                             &
           )
      MPH_ONFAILURE_RETURN(info)

      Call XMPH_PART_compute_interface_weight( mphs%lc_domain, info )
      MPH_ONFAILURE_RETURN(info)

      time_input2localsystem = MPI_Wtime() - time_input2localsystem

      !-------------------------------------------------------------------------
      ! [-] Save results & statistics
      !-------------------------------------------------------------------------

      Call MPH_mem_add2usage(mphs%mem,&
           XMPH_sm_sizeof(mphs%sls_domain%sm_A))
      mphs%rinfo (RINFO_TIMING_In2LcSys) = time_input2localsystem
      mphs%iinfo (IINFO_STRAT_PART     ) = mphs%ikeep( IKEEP_PART_Strategy )

    End Subroutine XMPH_distsystem_compute






  End Module XMPH_distsystem_mod
