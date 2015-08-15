! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"

!> testing module for dds_mod.
!!
!!
Module XMPH_test_dds_mod

  !* Modules *!

  Use XMPH_dds_mod
  Use test_utils_mod
  Use XMPH_test_fixture_mod

  !* No implicit typing *!
  Implicit None

  !* Access Specifiers *!
  Public :: XMPH_test_dds_init
  Public :: XMPH_test_dds_exit
  Public :: XMPH_test_dds_solve_sym
  Public :: XMPH_test_dds_solve_unsym
  Public :: XMPH_test_dds_solve_spd

  Private :: XMPH_test_dds_solve

Contains

  ! [+] routine : XMPH_test_dds_init ------------------------------------
  !
  !> Initialize the XMPH_test_sparse_matrix module
  !! 
  !! Namely, give default values to shared variables
  !!
  Subroutine XMPH_test_dds_init
    Call XMPH_test_fixture_init
  End Subroutine XMPH_test_dds_init


  ! [+] routine : XMPH_test_dds_exit ------------------------------------
  !
  !> Exit the XMPH_test_sparse_matrix module
  !! 
  !! Namely, free the memory used by the shared variables
  !!
  Subroutine XMPH_test_dds_exit
    Call XMPH_test_fixture_exit
  End Subroutine XMPH_test_dds_exit


  ! [+] routine : XMPH_test_dds_solve_unsym -----------------------------------------
  !
  !> testing the matrix vector product : A.x = y, with A unsymmetric 
  Subroutine XMPH_test_dds_solve_unsym

    Real, Parameter :: eps   = 1.e-6
    Call XMPH_test_dds_solve( XMPH_A1, XMPH_x1, XMPH_b1, eps )

  End Subroutine XMPH_test_dds_solve_unsym


  ! [+] routine : XMPH_test_dds_solve_sym -------------------------------------------
  !
  !> testing the matrix vector product : A.x = y, with A symmetric 
  Subroutine XMPH_test_dds_solve_sym
    Real, Parameter :: eps   = 1.e-6

    Call XMPH_test_dds_solve( XMPH_A2, XMPH_x2, XMPH_b2, eps )
  End Subroutine XMPH_test_dds_solve_sym


  ! [+] routine : XMPH_test_dds_solve_spd -------------------------------------------
  !
  !> testing the matrix vector product : A.x = y, with A spd
  Subroutine XMPH_test_dds_solve_spd
    Real, Parameter :: eps   = 1.e-6

    Call XMPH_test_dds_solve( XMPH_A3, XMPH_x3, XMPH_b3, eps )
  End Subroutine XMPH_test_dds_solve_spd



  ! [+] routine : XMPH_test_dds_solve ------------------------------------------------
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
  !!
  !!----
  !!
  !! @note  on output, it modify XMPH_test_result
  !! @warning  do not check the status of called subroutines
  !!
  Subroutine XMPH_test_dds_solve (A,x,b, eps)

    !* Arguments *!
    Type(XMPH_dense_matrix_t), Intent(in) :: A
    Type(XMPH_dense_matrix_t), Intent(in) :: x
    Type(XMPH_dense_matrix_t), Intent(in) :: b
    Real                , Intent(in) :: eps 

    !* Local variables *!
    Type(XMPH_dds_t) :: dds  ! dense direct solver
    Type(XMPH_dense_matrix_t       ) :: rhs  ! right hand side / solution
    Type(XMPH_dense_matrix_t       ) :: A_1  ! factors of A
    Integer :: status
    
    !-------------------------------------------------------------------------
    ! [1] Setup the fixture
    !-------------------------------------------------------------------------

    test_result = TEST_UNSET
    Call XMPH_dm_dup( A, A_1 , status )
    Call XMPH_dm_dup( b, rhs , status )

    !-------------------------------------------------------------------------
    ! [2] Perform the test
    !-------------------------------------------------------------------------
    
    Call XMPH_dds_init(dds,status)
    Call XMPH_dds_factorize( dds, A_1, status )
    Call XMPH_dds_solve( dds, A_1, rhs, status )

    ! Write (*,*) "rhs=", rhs%v
    Write(test_msg,*) " | y-b | =",&
         XMPH_norm2_diff(x%v, rhs%v, x%m), "(max=", eps,")"
    If ( eps <  XMPH_norm_inf(x%v, rhs%v, x%m ) )Then
       test_result = TEST_FAIL
    Else
       test_result = TEST_PASS
    End If

    !-------------------------------------------------------------------------
    ! [3] Tear Down
    !-------------------------------------------------------------------------

    Call XMPH_dds_exit(dds, status)
    Call XMPH_dm_free(rhs, status )

  End Subroutine XMPH_test_dds_solve


End Module XMPH_test_dds_mod
