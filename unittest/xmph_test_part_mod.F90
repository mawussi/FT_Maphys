! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"
!> testing module for XMPH_part_mod.
!!
!!
Module XMPH_test_part_mod

  !* Modules *!
  Use MPH_part_type
  Use XMPH_part_mod
  Use XMPH_sparse_matrix_mod
  Use mph_error_mod
  Use mph_dbg_mod
  Use test_utils_mod

  !* No implicit typing *!
  Implicit None

  !* Access Specifiers *!
  Public :: XMPH_test_part_init
  Public :: XMPH_test_part_exit
  Public :: XMPH_test_part_ordmat_metisnode
  Public :: XMPH_test_part_ordmat_metisedge
  Public :: XMPH_test_part_ordmat_metiswnode
  Public :: XMPH_test_part_ordmat_scotchnode
  Public :: XMPH_test_part_build_domains

  Private :: XMPH_test_part

  Private :: gen_poisson2D_5pts

  !* Private data
  Integer, Private, Parameter :: ORDER_MATRIX=1
  Integer, Private, Parameter :: BUILD_DOMAINS=2
  Integer, Private, Parameter :: DIST_MATRIX=3

  Integer, Private :: nd  ! number of domains
  Integer, Private, Parameter :: grid_size = 50
  Type(XMPH_sparse_matrix_t), Private, Save :: sm_Poisson ! 2D poisson 5 points matrix
  Character(len=MAPHYS_STRL), Private, Parameter :: Flname = &
       "XMPH_ARITHmphtest_part_mod.F90"

  Integer, Private :: mpirank
  Integer, Private :: mpisize

Contains
  Subroutine XMPH_test_part_init
    Use mph_log_mod
    Include "mpif.h"

    Call MPI_Comm_rank(MPI_COMM_WORLD,mpirank,test_result)
    ASSRT(test_result==MPI_SUCCESS)
    Call MPI_Comm_size(MPI_COMM_WORLD,mpisize,test_result)
    ASSRT(test_result==MPI_SUCCESS)
    
    If (mpirank == 0 ) Call gen_poisson2D_5pts
    !If (mpirank == 0 ) Call XMPH_sm_mmwrite(sm_Poisson,10,test_result)
    If (mpisize == 1 ) nd = 4
    If (mpisize /= 1 ) nd = mpisize

    Call MPH_dbg_init
    Call mph_log_init(6,6,6,MSG_ERROR,mpirank)

  End Subroutine XMPH_test_part_init

  Subroutine XMPH_test_part_exit
    Call XMPH_sm_free(sm_Poisson,test_result)
    ASSRT(test_result == 0)
    Call MPH_dbg_exit
  End Subroutine XMPH_test_part_exit


#ifdef HAVE_LIBMETIS
  !> test the metis node based partitioning strategy **
  Subroutine XMPH_test_part_ordmat_metisnode
    Call MPH_dbg_set_file("XMPH_ARITHpart_metisnode")
    Call XMPH_test_part(sm_Poisson,ORDER_MATRIX,nd,1,dbg_unit)
  End Subroutine XMPH_test_part_ordmat_metisnode

  !> test the metis edge based partitioning strategy 
  Subroutine XMPH_test_part_ordmat_metisedge
    Call MPH_dbg_set_file("XMPH_ARITHpart_metisedge")
    Call XMPH_test_part(sm_Poisson,ORDER_MATRIX,nd,2,dbg_unit)
  End Subroutine XMPH_test_part_ordmat_metisedge

  !> test the metis weighted node based partitioning strategy
  Subroutine XMPH_test_part_ordmat_metiswnode
    Call MPH_dbg_set_file("XMPH_ARITHpart_metiswnode")
    Call XMPH_test_part(sm_Poisson,ORDER_MATRIX,nd,3,dbg_unit)
  End Subroutine XMPH_test_part_ordmat_metiswnode

#else

  Subroutine XMPH_test_part_ordmat_metisnode
    test_result = TEST_SKIP
  End Subroutine XMPH_test_part_ordmat_metisnode

  Subroutine XMPH_test_part_ordmat_metisedge
    test_result = TEST_SKIP
  End Subroutine XMPH_test_part_ordmat_metisedge

  Subroutine XMPH_test_part_ordmat_metiswnode
    test_result = TEST_SKIP
  End Subroutine XMPH_test_part_ordmat_metiswnode

#endif

#ifdef HAVE_LIBSCOTCH
  !> test the scotch partitioning strategy
  Subroutine XMPH_test_part_ordmat_scotchnode
    Call MPH_dbg_set_file("XMPH_ARITHpart_scotchnode")
    Call XMPH_test_part(sm_Poisson,ORDER_MATRIX,nd,4,dbg_unit)
  End Subroutine XMPH_test_part_ordmat_scotchnode

  !> test the domains computation
  Subroutine XMPH_test_part_build_domains
    Call XMPH_test_part(sm_Poisson,BUILD_DOMAINS,nd,4,-1)
  End Subroutine XMPH_test_part_build_domains

  !> test the domains computation
  Subroutine XMPH_test_part_dist_matrix
    Call XMPH_test_part(sm_Poisson,DIST_MATRIX,nd,4,-1)
  End Subroutine XMPH_test_part_dist_matrix

#else
  Subroutine XMPH_test_part_ordmat_scotchnode
    test_result = TEST_SKIP
  End Subroutine XMPH_test_part_ordmat_scotchnode
  Subroutine XMPH_test_part_build_domains
    test_result = TEST_SKIP
  End Subroutine XMPH_test_part_build_domains
  Subroutine XMPH_test_part_dist_matrix
    test_result = TEST_SKIP
  End Subroutine XMPH_test_part_dist_matrix
#endif



  !> test if XMPH_PART routines  exit correctly
  !!
  !! the reordered matrix is written in "unit", if it is a valid file unit.
  !!
  !! @param [in] sm    the sparse matrix to reorder
  !! @param [in] nbdom the number of nested dissection to perform
  !! @param [in] strat the paritioning strategy
  !! @param [in] unit the sparse matrix reordered

  Subroutine XMPH_test_part(smin,job,nbdom,strat,unit)

    Include "mpif.h"
    Type(XMPH_sparse_matrix_t), Intent(in) :: smin
    Integer, Intent(in) :: job
    Integer, Intent(in) :: nbdom
    Integer, Intent(in) :: strat
    Integer, Intent(in) :: unit

    Integer :: status
    Integer :: min, max
    Integer :: comm

    Real(kind=8)                :: ana_timing(ANA_TIMING_SIZE)

    Type(XMPH_sparse_matrix_t)  :: sm
    Type(XMPH_sparse_matrix_t)  :: sm_local
    Type(maphys_matrix_graph_t) :: graph
    Type(maphys_binary_tree_t)  :: tree
    Type(maphys_domains_t) :: domains
    Type(maphys_domain_t) :: this_domain
    Type(maphys_rhs_partition_t) :: part_rhs

    ! Exit early
    If ((job >= DIST_MATRIX).And.(mpisize == 1 ))Then 
       ! avoid this test in sequential
       test_result = TEST_SKIP
       Write(test_msg,*) "Number of MPI Processes == 1"
       Return
    End If

    ! Init
    status = 0
    test_result = TEST_FAIL
    If (mpirank == 0 ) Call XMPH_sm_dup(sm,smin,status)

    ! Perform the jobs
    If (mpirank == 0)Then
       If (job >= ORDER_MATRIX) Then
          Call XMPH_PART_OrderGlobalMatrix   &      
               ( nbdom, strat, sm, ana_timing, &      
               & graph, tree, status )  
          If (status /= 0 ) Goto 9999
          If (unit >= 0) Call XMPH_sm_mmwrite(sm,unit,status)
       End If
       
       If (job >= BUILD_DOMAINS) Then
          Call XMPH_Part_build_domains &
               (sm, graph, tree, ana_timing,domains, status )
          If (status /= 0 ) Goto 9999
       End If
    End If

    If (job >= DIST_MATRIX) Then
       comm = MPI_COMM_WORLD
       Call XMPH_PART_DistGlobalMatrix  &
       (mpirank, comm , domains, &
       sm, graph, ana_timing, &
       part_rhs, this_domain, sm_local,      &
       status )
       If (status /= 0 ) Goto 9999
    End If

    If (status == 0 ) test_result = TEST_PASS

    ! Tear Down
9999 Continue

    ! error < 0 , warning > 0
    Call MPI_AllReduce(test_result,max,1,&
         MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,status)
    ASSRT(status == MPI_SUCCESS)
    test_result = max

    Call MPI_AllReduce(test_result,min,1,&
         MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,status)
    ASSRT(status == MPI_SUCCESS)
    If ( min < 0 ) test_result = min

    If (mpirank == 0 ) Call XMPH_sm_free(sm,status)

  End Subroutine XMPH_test_part


  !----------------------------------------------
  ! Generators
  !----------------------------------------------
  Subroutine gen_poisson2D_5pts
    Integer :: nnz, nD , nI, nIlow, nIupp
    Integer :: D_n , I_n, D_nnz, I_nnz
    Integer :: i, blk, idx

    nIlow = grid_size - 1
    nIupp = grid_size - 1
    nI = nIlow + nIupp
    I_n   = grid_size
    I_nnz = grid_size

    nD = grid_size
    D_n   = grid_size
    D_nnz = 3*grid_size -2

    nnz = nI * I_nnz + nD * D_nnz
    Call XMPH_sm_create(sm_Poisson, nnz, test_result)
    ASSRT(test_result == 0)

    sm_Poisson%m = grid_size * grid_size
    sm_Poisson%n = grid_size * grid_size

    ! all D blocs
    idx=0
    Do blk = 1, nD
       ! one D bloc

       ! diagonal
       Do i=1, D_n
          idx=idx+1
          sm_Poisson%i(idx) = i + (blk-1)*D_n
          sm_Poisson%j(idx) = i + (blk-1)*D_n
          sm_Poisson%v(idx) = 4.
       End Do

       ! upper
       Do i=1, D_n-1
          idx=idx+1 
          sm_Poisson%i(idx) = i   + (blk-1)*D_n
          sm_Poisson%j(idx) = i+1 + (blk-1)*D_n
          sm_Poisson%v(idx) = -1.
       End Do

       ! lower
       Do i=2, D_n
          idx=idx+1
          sm_Poisson%i(idx) = i    + (blk-1)*D_n
          sm_Poisson%j(idx) = i-1  + (blk-1)*D_n
          sm_Poisson%v(idx) = -1.
       End Do
    End Do

    ASSRT( idx == nD*D_nnz )

    ! All uppers I blocs
    Do blk = 1, nIupp 

       ! diagonal
       Do i=1, I_n
          idx=idx+1
          sm_Poisson%i(idx) = i + (blk-1)*I_n
          sm_Poisson%j(idx) = i + (blk  )*I_n
          sm_Poisson%v(idx) = -1.
       End Do

    End Do
    
    ! All lower I blocs
    Do blk = 2, nIlow + 1

       ! diagonal
       Do i=1, I_n
          idx=idx+1
          sm_Poisson%i(idx) = i + (blk-1)*I_n
          sm_Poisson%j(idx) = i + (blk-2)*I_n
          sm_Poisson%v(idx) = -1.
       End Do

    End Do

    ASSRT( idx == sm_Poisson%nnz )

  End Subroutine gen_poisson2D_5pts



End Module XMPH_test_part_mod
