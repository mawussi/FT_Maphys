! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"

!> testing sds_mod.
Module XMPH_test_sds_mod

  !* Modules *!

  Use XMPH_sds_mod
  Use XMPH_dense_matrix_mod
  Use test_utils_mod
  Use XMPH_test_fixture_mod

  !* No implicit typing *!
  Implicit None

  !* Local constants *!
  Character(len=MAPHYS_STRL), Private, Parameter :: FLNAME= &
       "XMPH_ARITHmph_test_sds_mod.F90"


  !* Access Specifiers *!
  Public :: XMPH_test_sds_init
  Public :: XMPH_test_sds_exit
  Public :: XMPH_test_sds_mumps_solve_sym
  Public :: XMPH_test_sds_mumps_solve_unsym
  Public :: XMPH_test_sds_mumps_solve_spd
  Public :: XMPH_test_sds_pastix_solve_sym
  Public :: XMPH_test_sds_pastix_solve_unsym
  Public :: XMPH_test_sds_pastix_solve_spd

  Public :: XMPH_test_sds_mumps_schur_unsym
  Public :: XMPH_test_sds_mumps_schur_sym
  Public :: XMPH_test_sds_mumps_schur_spd
  Public :: XMPH_test_sds_pastix_schur_unsym
  Public :: XMPH_test_sds_pastix_schur_sym
  Public :: XMPH_test_sds_pastix_schur_spd

  Private :: XMPH_test_sds_solve
  Private :: XMPH_test_sds_schur
  
Contains


  ! [+] routine : XMPH_test_sds_init ------------------------------------
  !
  !> Initialize the XMPH_test_sparse_matrix module
  !! 
  !! Namely, give default values to shared variables
  !!
  Subroutine XMPH_test_sds_init
    Call XMPH_test_fixture_init
  End Subroutine XMPH_test_sds_init


  ! [+] routine : XMPH_test_sds_exit ------------------------------------
  !
  !> Exit the XMPH_test_sparse_matrix module
  !! 
  !! Namely, free the memory used by the shared variables
  !!
  Subroutine XMPH_test_sds_exit
    Call XMPH_test_fixture_exit
  End Subroutine XMPH_test_sds_exit



  ! [+] routine : XMPH_test_sds_mumps_solve_unsym ------------------------------
  !
  !> testing the matrix vector product : A.x = y, with A unsymmetric 
  Subroutine XMPH_test_sds_mumps_solve_unsym
    Call XMPH_test_sds_solve( XMPH_spA1, XMPH_x1, XMPH_b1, 1.e-6, SDS_IsMUMPS )
  End Subroutine XMPH_test_sds_mumps_solve_unsym


  ! [+] routine : XMPH_test_sds_mumps_solve_sym --------------------------------
  !
  !> testing the matrix vector product : A.x = y, with A symmetric 
  Subroutine XMPH_test_sds_mumps_solve_sym
    Call XMPH_test_sds_solve( XMPH_spA2, XMPH_x2, XMPH_b2, 1.e-6, SDS_IsMUMPS )
  End Subroutine XMPH_test_sds_mumps_solve_sym


  ! [+] routine : XMPH_test_sds_mumps_solve_spd --------------------------------
  !
  !> testing the matrix vector product : A.x = y, with A spd
  Subroutine XMPH_test_sds_mumps_solve_spd
    Call XMPH_test_sds_solve( XMPH_spA3, XMPH_x3, XMPH_b3, 1.e-6, SDS_IsMUMPS )
  End Subroutine XMPH_test_sds_mumps_solve_spd


  ! [+] routine : XMPH_test_sds_pastix_solve_unsym -----------------------------
  !
  !> testing the matrix vector product : A.x = y, with A unsymmetric 
  Subroutine XMPH_test_sds_pastix_solve_unsym

    Call XMPH_test_sds_solve( XMPH_spA1, XMPH_x1, XMPH_b1, 1.e-6, SDS_IsPASTIX )

  End Subroutine XMPH_test_sds_pastix_solve_unsym


  ! [+] routine : XMPH_test_sds_pastix_solve_sym -------------------------------
  !
  !> testing the matrix vector product : A.x = y, with A symmetric 
  Subroutine XMPH_test_sds_pastix_solve_sym
    Call XMPH_test_sds_solve( XMPH_spA2, XMPH_x2, XMPH_b2, 1.e-6, SDS_IsPASTIX )
  End Subroutine XMPH_test_sds_pastix_solve_sym

  ! [+] routine : XMPH_test_sds_pastix_solve_spd -------------------------------
  !
  !> testing the matrix vector product : A.x = y, with A spd
  Subroutine XMPH_test_sds_pastix_solve_spd
    Call XMPH_test_sds_solve( XMPH_spA3, XMPH_x3, XMPH_b3, 1.e-6, SDS_IsPASTIX )
  End Subroutine XMPH_test_sds_pastix_solve_spd

  ! [+] routine : XMPH_test_sds_mumps_schur_unsym ------------------------------
  !
  !> testing the schur complement computation with A an unsym matrix
  Subroutine XMPH_test_sds_mumps_schur_unsym
    Call XMPH_test_sds_schur &
         ( XMPH_spA1, XMPH_nS1, XMPH_S1list, XMPH_S1, 1.e-6, SDS_IsMUMPS )
  End Subroutine XMPH_test_sds_mumps_schur_unsym

  ! [+] routine : XMPH_test_sds_mumps_schur_sym --------------------------------
  !
  !> testing the schur complement computation with A an sym matrix
  Subroutine XMPH_test_sds_mumps_schur_sym
    Call XMPH_test_sds_schur &
         ( XMPH_spA2, XMPH_nS2, XMPH_S2list, XMPH_S2, 1.e-6, SDS_IsMUMPS )
  End Subroutine XMPH_test_sds_mumps_schur_sym

  ! [+] routine : XMPH_test_sds_mumps_schur_spd --------------------------------
  !
  !> testing the schur complement computation with A an spd matrix
  Subroutine XMPH_test_sds_mumps_schur_spd
    Call XMPH_test_sds_schur &
         ( XMPH_spA3, XMPH_nS3, XMPH_S3list, XMPH_S3, 1.e-6, SDS_IsMUMPS )
  End Subroutine XMPH_test_sds_mumps_schur_spd


  ! [+] routine : XMPH_test_sds_pastix_schur_unsym -----------------------------
  !
  !> testing the schur complement computation with A an unsym matrix
  Subroutine XMPH_test_sds_pastix_schur_unsym
    Call XMPH_test_sds_schur &
         ( XMPH_spA1, XMPH_nS1, XMPH_S1list, XMPH_S1, 1.e-6, SDS_IsPASTIX )
  End Subroutine XMPH_test_sds_pastix_schur_unsym

  ! [+] routine : XMPH_test_sds_pastix_schur_sym -------------------------------
  !
  !> testing the schur complement computation with A an sym matrix
  Subroutine XMPH_test_sds_pastix_schur_sym
    Call XMPH_test_sds_schur &
         ( XMPH_spA2, XMPH_nS2, XMPH_S2list, XMPH_S2, 1.e-6, SDS_IsPASTIX )
  End Subroutine XMPH_test_sds_pastix_schur_sym

  ! [+] routine : XMPH_test_sds_pastix_schur_spd -------------------------------
  !
  !> testing the schur complement computation with A an spd matrix
  Subroutine XMPH_test_sds_pastix_schur_spd
    Call XMPH_test_sds_schur &
         ( XMPH_spA3, XMPH_nS3, XMPH_S3list, XMPH_S3, 1.e-6, SDS_IsPASTIX )
  End Subroutine XMPH_test_sds_pastix_schur_spd

  ! [+] routine : XMPH_test_sds_solve ------------------------------------------
  !
  !> test that dds can solve the linear system "A.x = b"
  !!
  !! Check that ( factorize(A) . b ~= x ) 
  !!
  !!----
  !!
  !! @param [in]  A   the matrix of the linear system
  !! @param [in]  x   the solution of the linear system
  !! @param [in]  b   the right hand side
  !! @param [in]  eps the acceptable error on | x_computed - x | (with norm infinite)
  !! @param [in]  slv the selected solver (2= MUMPS 3=PASTIX)
  !!
  !!----
  !!
  !! @note  on output, it modify test_result
  !! @warning  do not check the status of called subroutines
  !!
  Subroutine XMPH_test_sds_solve (spA,x,b, eps, slv )

    !* Module(s) & co. *!
    Use XMPH_sparse_matrix_mod
    Implicit None
    Include "mpif.h"

    !* Arguments *!
    Type(XMPH_sparse_matrix_t), Intent(in) :: spA
    Type(XMPH_dense_matrix_t), Intent(in) :: x
    Type(XMPH_dense_matrix_t), Intent(in) :: b
    Real                , Intent(in) :: eps 
    Integer             , Intent(in) :: slv

    !* Local variables *!
    Type(XMPH_sds_t) :: sds  ! dense direct solver
    Type(XMPH_dense_matrix_t) :: rhs  ! right hand side / solution
    Real(kind=8) :: cmp
    Integer :: status
    
    !-------------------------------------------------------------------------
    ! [1] Setup the fixture
    !-------------------------------------------------------------------------

    test_result = TEST_UNSET
    Call XMPH_sds_select(slv,sds,status)
    If( status /= 0 )Then
       test_result = TEST_SKIP
       Return
    End If
    Call XMPH_dm_dup( b, rhs , status )

    !-------------------------------------------------------------------------
    ! [2] Perform the test
    !-------------------------------------------------------------------------
    
    Call XMPH_sds_set_MPIcommunicator(sds,MPI_COMM_SELF,status)
    Call XMPH_sds_set_matrix(sds,spA,status)

    Call XMPH_sds_analyze(sds,status)
    Call XMPH_sds_factorize(sds,status)
    Call XMPH_sds_solve_RHS(sds,rhs%v,rhs%n,rhs%ld,status )
    Call XMPH_sds_finalize(sds, status )

    ! Write (*,*) "rhs=", rhs%v
    If (status == 0 ) Then
       cmp = XMPH_norm2_diff(x%v, rhs%v, x%m)
    Else
       cmp = 9999.d0
    End If

    Write(test_msg,*) " | y-b | =",cmp , "(max=", eps,")"
    If ( eps <  cmp )Then
       test_result = TEST_FAIL
    Else
       test_result = TEST_PASS
    End If

    !-------------------------------------------------------------------------
    ! [3] Tear Down
    !-------------------------------------------------------------------------

    Call XMPH_dm_free(rhs, status )

  End Subroutine XMPH_test_sds_solve


  ! [+] routine : XMPH_test_sds_schur ------------------------------------------
  !
  !> test if sds can compute the schur complement "A.x = b"
  !!
  !! Check that ( S - Stheoric ~= 0 ) 
  !!
  !!----
  !!
  !! @param [in]  spA       the matrix of the linear system
  !! @param [in]  nS        the order of the schur complement
  !! @param [in]  Slist     the indices of in spA of the schur complement
  !! @param [in]  Stheoric  the theoric schur complement
  !! @param [in]  eps the acceptable error on | x_computed - x | (with norm infinite)
  !! @param [in]  slv the selected solver (2= MUMPS 3=PASTIX)
  !!
  !!----
  !!
  !! @note  on output, it modify test_result
  !! @warning  do not check the status of called subroutines
  !!
  Subroutine XMPH_test_sds_schur (spA,nS, Slist , Stheoric, eps, slv )

    !* Module(s) & co. *!
    Use XMPH_sparse_matrix_mod
    Use mph_error_mod
    Implicit None
    Include "mpif.h"

    !* Arguments *!
    Type(XMPH_sparse_matrix_t), Intent(in) :: spA
    Integer             , Intent(in) :: nS
    Integer             , Intent(in) :: Slist (:)
    Type(XMPH_dense_matrix_t), Intent(in) :: Stheoric
    Real                , Intent(in) :: eps 
    Integer             , Intent(in) :: slv

    !* Local variables *!
    Type(XMPH_sds_t) :: sds  ! dense direct solver
    Type(XMPH_dense_matrix_t)         :: S    ! the computed schur
    Real(kind=8) :: cmp
    Integer :: status
    
    !-------------------------------------------------------------------------
    ! [1] Setup the fixture
    !-------------------------------------------------------------------------

    test_result = TEST_UNSET
    Call XMPH_sds_select(slv,sds,status)
    If( status /= 0 )Then
       test_result = TEST_SKIP
       Return
    End If

    !-------------------------------------------------------------------------
    ! [2] Perform the test
    !-------------------------------------------------------------------------
  
    Call XMPH_sds_set_MPIcommunicator(sds,MPI_COMM_SELF,status)
    Call XMPH_sds_set_matrix(sds,spA,status)

    Call XMPH_sds_set_schurlist(sds,nS,Slist,status)
    ASSRT(status==0)
    Call XMPH_sds_analyze(sds,status)
    ASSRT(status==0)
    Call XMPH_sds_factorize(sds,status)
    ASSRT(status==0)
    Call XMPH_sds_get_schur(sds,S,status)
    ASSRT(status==0)
    Call XMPH_dm_unsymmetrize(S) 
    
    ! Write(6,*) ""
    ! Write(6,*) "S="
    ! Call XMPH_dm_mmwrite(S,6,status)
    ! Write(6,*) "Sref="
    ! Call XMPH_dm_mmwrite(Stheoric,6,status)
    If (status == 0 ) Then
       cmp = XMPH_norm2_diff(S%v, Stheoric%v, nS*nS)
    Else
       cmp = 9999.d0
    End If

    Write(test_msg,*) " | S-Sref | =", cmp , "(max=", eps,")"
    If ( eps <  cmp )Then
       test_result = TEST_FAIL
    Else
       test_result = TEST_PASS
    End If

    !-------------------------------------------------------------------------
    ! [3] Tear Down
    !-------------------------------------------------------------------------

    Call XMPH_sds_finalize(sds, status )
    Call XMPH_dm_nullify(S, status )

  End Subroutine XMPH_test_sds_schur



End Module XMPH_test_sds_mod

