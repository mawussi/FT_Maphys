! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"

!> testing ilu_mod.
Module XMPH_test_ilu_mod

  !* Modules *!

  ! Use XMPH_ilu_mod
  Use XMPH_pilut_mod
  Use XMPH_dense_matrix_mod
  Use test_utils_mod
  Use XMPH_test_fixture_mod

  !* No implicit typing *!
  Implicit None

  !* Local constants *!
  Character(len=MAPHYS_STRL), Private, Parameter :: FLNAME= &
       "XMPH_ARITHmph_test_ilu_mod.F90"


  !* Access Specifiers *!
  Public :: XMPH_test_ilu_init
  Public :: XMPH_test_ilu_exit

  Public :: XMPH_test_pilut_schur_unsym
  ! Public :: XMPH_test_pilut_schur_sym
  ! Public :: XMPH_test_pilut_schur_spd

  Private :: XMPH_test_ilu_schur
  
Contains


  ! [+] routine : XMPH_test_ilu_init ------------------------------------
  !
  !> Initialize the XMPH_test_sparse_matrix module
  !! 
  !! Namely, give default values to shared variables
  !!
  Subroutine XMPH_test_ilu_init
    Call XMPH_test_fixture_init
  End Subroutine XMPH_test_ilu_init


  ! [+] routine : XMPH_test_ilu_exit ------------------------------------
  !
  !> Exit the XMPH_test_sparse_matrix module
  !! 
  !! Namely, free the memory used by the shared variables
  !!
  Subroutine XMPH_test_ilu_exit
    Call XMPH_test_fixture_exit
  End Subroutine XMPH_test_ilu_exit

  ! [+] routine : XMPH_test_pilut_schur_unsym ------------------------------
  !
  !> testing the schur complement computation with A an unsym matrix
  Subroutine XMPH_test_pilut_schur_unsym

    MPH_INT                :: ilims (2)
    Real(KIND=XMPH_FLOATKIND) :: rlims (2)

    ilims = (/ 100, 100 /) 
    rlims = (/ 10._XMPH_FLOATKIND, 10._XMPH_FLOATKIND /) 

    Call XMPH_test_ilu_schur &
         ( XMPH_spA1, XMPH_nS1, ilims, rlims, XMPH_SpS1, 1.e-6 )

  End Subroutine XMPH_test_pilut_schur_unsym

  ! [+] routine : XMPH_test_pilut_schur_sym --------------------------------
  !
  !> testing the schur complement computation with A an sym matrix
  ! Subroutine XMPH_test_pilut_schur_sym
  !   Call XMPH_test_ilu_schur &
  !        ( XMPH_spA2, XMPH_nS2, XMPH_S2list, XMPH_S2, 1.e-6, ILU_IsMUMPS )
  ! End Subroutine XMPH_test_pilut_schur_sym

  ! [+] routine : XMPH_test_pilut_schur_spd --------------------------------
  !
  !> testing the schur complement computation with A an spd matrix
  ! Subroutine XMPH_test_pilut_schur_spd
  !   Call XMPH_test_ilu_schur &
  !        ( XMPH_spA3, XMPH_nS3, XMPH_S3list, XMPH_S3, 1.e-6, ILU_IsMUMPS )
  ! End Subroutine XMPH_test_pilut_schur_spd

  ! [+] routine : XMPH_test_ilu_schur ------------------------------------------
  !
  !> test if ilu can compute the schur complement "A.x = b"
  !!
  !! Check that ( S - Stheoric ~= 0 ) 
  !!
  !!----
  !!
  !! @param [in]  spA       the matrix of the linear system
  !! @param [in]  nS        the order of the schur complement
  !! @param [in]  ilims     the limits on sizes 
  !!       ilims(1) : maximum number of LU entries
  !!       ilims(2) : maximum number of schur entries
  !! @param [in]  rlims     the limits on values
  !!       ilims(1) : tolerance on LU factors
  !!       ilims(2) : tolerance on Schur complement.
  !! @param [in]  Stheoric  the theoric schur complement
  !! @param [in]  eps the acceptable error on | x_computed - x | (with norm infinite)
  !!
  !!----
  !!
  !! @note  on output, it modify test_result
  !! @warning  do not check the status of called subroutines
  !!
  Subroutine XMPH_test_ilu_schur (spA,nS, ilims, rlims , Stheoric, eps )

    !* Module(s) & co. *!
    Use XMPH_sparse_matrix_mod
    Use mph_error_mod
    Use XMPH_pilut_mod
    Implicit None
    Include "mpif.h"

    !* Arguments *!
    Type(XMPH_sparse_matrix_t) , Intent(in) :: spA
    Integer                    , Intent(in) :: nS
    MPH_INT                 , Intent(in) :: ilims (2)
    Real(KIND=XMPH_FLOATKIND)  , Intent(in) :: rlims (2)
    Type(XMPH_sparse_matrix_t) , Intent(in) :: Stheoric
    Real                       , Intent(in) :: eps 

    !* Local variables *!
    Type(XMPH_pilut_t)         :: this ! the ilu solver
    Type(XMPH_sparse_matrix_t) :: S    ! the computed schur
    Real(kind=8)               :: cmp
    Integer                    :: status
    
    !-------------------------------------------------------------------------
    ! [1] Setup the fixture
    !-------------------------------------------------------------------------

    test_result = TEST_UNSET
    Call XMPH_pilut_init(this)

    Call XMPH_pilut_set_matrix(this,spA,status)
    If (status /= 0) test_result = TEST_FAIL
    If (test_result == TEST_FAIL ) Goto 9999

    Call XMPH_pilut_set_ILULimits(this,&
         ilims(1),rlims(1),&
         ilims(2),rlims(2),status)
    If (status /= 0) test_result = TEST_FAIL
    If (test_result == TEST_FAIL ) Goto 9999

    Call XMPH_pilut_set_schurOrder(this,nS,status)
    If (status /= 0) test_result = TEST_FAIL
    If (test_result == TEST_FAIL ) Goto 9999

    !-------------------------------------------------------------------------
    ! [2] Perform the test
    !-------------------------------------------------------------------------

    Call XMPH_pilut_factorize(this,status)
    If (status /= 0) test_result = TEST_FAIL
    If (test_result == TEST_FAIL ) Goto 9999

    Call XMPH_pilut_get_schur(this,S,status)
    If (status /= 0) test_result = TEST_FAIL
    If (test_result == TEST_FAIL ) Goto 9999

    test_result = TEST_PASS
    Write(test_msg,*) "check based on return code only"

    !-------------------------------------------------------------------------
    ! [3] Tear down
    !-------------------------------------------------------------------------

9999 Continue    

    Call XMPH_pilut_exit(this)
    If (status == 0 ) Call XMPH_sm_free(S,status)

  End Subroutine XMPH_test_ilu_schur



End Module XMPH_test_ilu_mod

