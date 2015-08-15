! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
!> testing module for the module dense_matrix_mod.
Module XMPH_test_dense_matrix_mod

  !* Modules *!
  Use test_utils_mod
  Use XMPH_dense_matrix_mod 
  Use XMPH_test_fixture_mod
  Implicit None

  !* Access specifiers for routines *!
  Public :: XMPH_test_dm_init
  Public :: XMPH_test_dm_exit
  Public :: XMPH_test_dm_matvect_unsym
  Public :: XMPH_test_dm_matvect_sym
  Public :: XMPH_test_dm_matvect_spd

  Private :: XMPH_test_dm_matvect

Contains

  ! [+] routine : XMPH_test_dm_init ------------------------------------
  !
  !> Initialize the XMPH_test_dm module
  !! 
  !! Namely, give default values to shared variables
  !!
  Subroutine XMPH_test_dm_init
    Call XMPH_test_fixture_init
  End Subroutine XMPH_test_dm_init


  ! [+] routine : XMPH_test_dm_exit ------------------------------------
  !
  !> Exit the XMPH_test_dm module
  !! 
  !! Namely, free the memory used by the shared variables
  !!
  Subroutine XMPH_test_dm_exit
    Call XMPH_test_fixture_exit
  End Subroutine XMPH_test_dm_exit


  ! [+] routine : XMPH_test_dm_matvect_unsym ------------------------------------------
  !
  !> testing the matrix vector product : A.x = y, with A unsymmetric 
  Subroutine XMPH_test_dm_matvect_unsym
    
    Call XMPH_test_dm_matvect (XMPH_A1,XMPH_x1,XMPH_b1,1.e-6)

  End Subroutine XMPH_test_dm_matvect_unsym


  ! [+] routine : XMPH_test_dm_matvect_sym -------------------------------------------
  !
  !> testing the matrix vector product : A.x = y, with A symmetric 
  Subroutine XMPH_test_dm_matvect_sym

    Call XMPH_test_dm_matvect (XMPH_A2,XMPH_x2,XMPH_b2,1.e-6)

  End Subroutine XMPH_test_dm_matvect_sym


  ! [+] routine : XMPH_test_dm_matvect_sym -------------------------------------------
  !
  !> testing the matrix vector product : A.x = y, with A symmetric 
  Subroutine XMPH_test_dm_matvect_spd

    Call XMPH_test_dm_matvect (XMPH_A3,XMPH_x3,XMPH_b3,1.e-6)

  End Subroutine XMPH_test_dm_matvect_spd


  ! [+] routine : XMPH_test_dm_matvect --------------------------------------------
  !
  !> test that dense matrix product "A.x = b"
  !!
  !! Check that ( A . x ~= b ) 
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
  !! @note  on output, it modify test_result
  !! @warning  do not check the status of called subroutines
  !!
  Subroutine XMPH_test_dm_matvect( A, x, b, eps )

    !* Arguments *!
    Type(XMPH_dense_matrix_t), Intent(in) :: A
    Type(XMPH_dense_matrix_t), Intent(in) :: x
    Type(XMPH_dense_matrix_t), Intent(in) :: b
    Real                , Intent(in) :: eps 

    !* Local variables *!
    Type(XMPH_dense_matrix_t       ) :: y  ! computed right hand side 
    Integer :: status
    
    !-------------------------------------------------------------------------
    ! [1] Setup the fixture
    !-------------------------------------------------------------------------

    test_result = TEST_UNSET

    Call XMPH_dm_create(y, x%m, x%n, x%ld , status )
    If ( status /= 0 ) Return

    !-------------------------------------------------------------------------
    ! [2] Perform the test
    !-------------------------------------------------------------------------

    Call XMPH_dm_vectorproduct(A, x, y, status )
    If ( status /= 0 ) Return

    !Write (*,*) "y=", y%v     
    Write(test_msg,*) " | y-b | =", &
         XMPH_norm2_diff(b%v, y%v, b%m), "(max=", eps,")"

    If ( eps <  XMPH_norm_inf(b%v, y%v, b%m ) )Then
       test_result = TEST_FAIL
    Else
       test_result = TEST_PASS
    End If

    !-------------------------------------------------------------------------
    ! [3] Tear Down
    !-------------------------------------------------------------------------

    Call XMPH_dm_free(y, status )
    If ( status /= 0 ) Return

  End Subroutine XMPH_test_dm_matvect

End Module XMPH_test_dense_matrix_mod
