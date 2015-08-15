! Warning: XMPH_GENFILE_COMMENT
  !> test the modules
  !! @note we only test here the interfaces.
  !! @author Yohan Lee-tin-yien
  Program XMPH_test_all

    Use test_utils_mod
    Use XMPH_test_dense_matrix_mod
    Use XMPH_test_sparse_matrix_mod
    Use XMPH_test_dds_mod
    Use XMPH_test_sds_mod
    Use XMPH_test_ilu_mod
    Use XMPH_test_part_mod

    Implicit None
    Include "mpif.h"
    Integer :: ierror
    Integer :: rank

    ! Init MPI
    Call MPI_Init(ierror)
    Call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierror )
    
    If ( rank == 0 ) Then

       ! dense_matrix
       Call XMPH_test_dm_init
       Call Test( XMPH_test_dm_matvect_unsym  , "XMPH_dm.matvect() - unsym")
       Call Test( XMPH_test_dm_matvect_sym    , "XMPH_dm.matvect() -   sym")
       Call Test( XMPH_test_dm_matvect_spd    , "XMPH_dm.matvect() -   spd")
       Call XMPH_test_dm_exit

       ! sparse_matrix
       Call XMPH_test_sm_init
       Call Test( XMPH_test_sm_ind2ptr       , "XMPH_sm.ind2ptr() ")
       Call Test( XMPH_test_sm_ptr2ind       , "XMPH_sm.ptr2ind() ")
       Call Test( XMPH_test_sm_qsort         , "XMPH_sm.qsort() ")
       Call Test( XMPH_test_sm_assemble      , "XMPH_sm.assemble() ")
       Call Test( XMPH_test_sm_symStruct     , "XMPH_sm.symStruct() ")
       Call Test( XMPH_test_sm_matvect_unsym , "XMPH_sm.matvect() - unsym")
       Call Test( XMPH_test_sm_matvect_sym   , "XMPH_sm.matvect() -   sym")
       Call Test( XMPH_test_sm_matvect_spd   , "XMPH_sm.matvect() -   spd")
       Call XMPH_test_sm_exit

       ! dense direct solver (LAPACK)
       Call XMPH_test_dds_init
       Call Test( XMPH_test_dds_solve_unsym  , "XMPH_dds.solve() - unsym")
       Call Test( XMPH_test_dds_solve_sym    , "XMPH_dds.solve() -   sym")
       Call Test( XMPH_test_dds_solve_spd    , "XMPH_dds.solve() -   spd")
       Call XMPH_test_dds_exit


       ! sparse direct solver (MUMPS/PASTIX)
       Call XMPH_test_sds_init
       Call Test( XMPH_test_sds_mumps_solve_unsym  , "sds_mumps.solve() - unsym")
       Call Test( XMPH_test_sds_mumps_solve_sym    , "sds_mumps.solve() -   sym")
       Call Test( XMPH_test_sds_mumps_solve_spd    , "sds_mumps.solve() -   spd")
       Call Test( XMPH_test_sds_pastix_solve_unsym , "sds_pastix.solve()- unsym")
       Call Test( XMPH_test_sds_pastix_solve_sym   , "sds_pastix.solve()-   sym")
       Call Test( XMPH_test_sds_pastix_solve_spd   , "sds_pastix.solve()-   spd")
       Call Test( XMPH_test_sds_mumps_schur_unsym  , "sds_mumps.schur() - unsym")
       Call Test( XMPH_test_sds_mumps_schur_sym    , "sds_mumps.schur() -   sym")
       Call Test( XMPH_test_sds_mumps_schur_spd    , "sds_mumps.schur() -   spd")
       Call Test( XMPH_test_sds_pastix_schur_unsym , "sds_pastix.schur()- unsym")
       Call Test( XMPH_test_sds_pastix_schur_sym   , "sds_pastix.schur()-   sym")
       Call Test( XMPH_test_sds_pastix_schur_spd   , "sds_pastix.schur()-   spd")
       Call XMPH_test_sds_exit

       ! ilu solver (pilut)
       Call XMPH_test_ilu_init
       Call Test( XMPH_test_pilut_schur_unsym  , "pilut.schur() - unsym")
       Call XMPH_test_ilu_exit
    End If
  
    ! Parallel tests
    ! Partitioner
    Call XMPH_test_part_init
    Call MPI_Barrier(MPI_COMM_WORLD,ierror)
    Call Test(XMPH_test_part_ordmat_metisnode,  "part.order_matrix()-  metis Node")
    Call MPI_Barrier(MPI_COMM_WORLD,ierror)
    Call Test(XMPH_test_part_ordmat_metisedge,  "part.order_matrix()-  metis Edge")
    Call MPI_Barrier(MPI_COMM_WORLD,ierror)
    Call Test(XMPH_test_part_ordmat_metiswnode, "part.order_matrix()- metis WNode")
    Call MPI_Barrier(MPI_COMM_WORLD,ierror)
    Call Test(XMPH_test_part_ordmat_scotchnode, "part.order_matrix()- scotch Node")
    Call MPI_Barrier(MPI_COMM_WORLD,ierror)
    Call Test(XMPH_test_part_build_domains,    "part.build_domains()- scotch Node")
    Call MPI_Barrier(MPI_COMM_WORLD,ierror)
    Call Test(XMPH_test_part_dist_matrix,       "part.build_matrix()- scotch Node")
    Call MPI_Barrier(MPI_COMM_WORLD,ierror)
    Call XMPH_test_part_exit

    ! Print the summary of all tests.
    If ( rank == 0 ) Call test_summarize

    ! End MPI
    Call MPI_Finalize(ierror)

  End Program XMPH_test_all


