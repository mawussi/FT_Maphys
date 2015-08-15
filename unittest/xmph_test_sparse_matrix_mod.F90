! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"

!> testing module for the module sparse_matrix_mod.
!!
!!
Module XMPH_test_sparse_matrix_mod
  Use XMPH_sparse_matrix_mod 
  Use test_utils_mod
  Use XMPH_test_fixture_mod
  Implicit None

  !* Private variables
  Integer, Parameter, Private :: ref_indsize = 3
  Integer, Parameter, Private :: ref_ptrsize = 3
  Integer, Parameter, Private :: ref_ind(ref_indsize) = (/ 1,1,2 /)
  Integer, Parameter, Private :: ref_ptr(ref_ptrsize) = (/ 1,3,4 /)

  !* Access Specifiers to routine *!
  Public :: XMPH_test_sm_init
  Public :: XMPH_test_sm_exit

  Public :: XMPH_test_sm_qsort
  Public :: XMPH_test_sm_assemble
  Public :: XMPH_test_sm_symStruct

  Public :: XMPH_test_sm_matvect_unsym
  Public :: XMPH_test_sm_matvect_sym
  Public :: XMPH_test_sm_matvect_spd
  Private :: XMPH_test_sm_matvect


Contains

  ! [+] routine : XMPH_test_sm_init ------------------------------------
  !
  !> Initialize the XMPH_test_sm module
  !! 
  !! Namely, give default values to shared variables
  !!
  Subroutine XMPH_test_sm_init
    Call XMPH_test_fixture_init
  End Subroutine XMPH_test_sm_init


  ! [+] routine : XMPH_test_sm_exit ------------------------------------
  !
  !> Exit the XMPH_test_sm module
  !! 
  !! Namely, free the memory used by the shared variables
  !!
  Subroutine XMPH_test_sm_exit
    Call XMPH_test_fixture_exit
  End Subroutine XMPH_test_sm_exit


  ! [+] routine : XMPH_test_sm_ind2ptr --------------------------------------
  !
  !> testing getting the pointer array from a list of indices
  Subroutine XMPH_test_sm_ind2ptr

    Integer :: ptr(ref_ptrsize)
    Integer :: maxdiff

    ! Setup
    test_result = TEST_UNSET

    ! Run
    Call XMPH_sm_ind2ptr(ref_indsize,ref_ind, ref_ptrsize, ptr )

    ! Check 
    ptr = ptr-ref_ptr
    maxdiff = MaxVal( ptr )
    If (maxdiff == 0 ) test_result = TEST_PASS
    If (maxdiff /= 0 ) test_result = TEST_FAIL

    Write(test_msg,*) "maxval(res-ref)=", maxdiff

  End Subroutine XMPH_test_sm_ind2ptr


  ! [+] routine : XMPH_test_sm_ptr2ind --------------------------------------
  !
  !> testing getting the list of indices from a pointer array
  Subroutine XMPH_test_sm_ptr2ind

    Integer :: ind(ref_indsize)
    Integer :: maxdiff

    ! Setup
    test_result = TEST_UNSET

    ! Run
    Call XMPH_sm_ptr2ind(ref_ptrsize,ref_ptr, ind )

    ! Check 
    ind = ind-ref_ind
    maxdiff = MaxVal( ind )
    If (maxdiff == 0 ) test_result = TEST_PASS
    If (maxdiff /= 0 ) test_result = TEST_FAIL

    Write(test_msg,*) "maxval(res-ref)=", maxdiff

  End Subroutine XMPH_test_sm_ptr2ind


  ! [+] routine : test_maphys_qsort ----------------------------------------
  !
  !> testing sorting with qsort an array
  Subroutine XMPH_test_sm_qsort

    MPH_INT ::  ref(6) = (/ 1, 2, 3, 4, 5, 6 /)
    MPH_INT ::  permref(6) = (/ 1, 5, 3, 4, 6, 2 /)
    MPH_INT, Pointer ::    i(:) 
    MPH_INT, Pointer :: perm(:) 
    Integer    :: maxdiff
    Integer    :: info

    ! Setup
    test_result = TEST_UNSET
    Allocate(i(6))
    Allocate(perm(6))
    i = (/ 1, 6, 3, 4, 2, 5 /)
    
    ! Run
    Call XMPH_qsort(6, i, perm, info)

    ! Check

    If (info >= 0 )  test_result = TEST_PASS
    If (info <  0 )  test_result = TEST_FAIL
    test_msg = "routine status checked"

    If ( test_result == TEST_PASS )Then
       ! print *, "i   = ", i
       ! print *, "ref = ", ref
       i = i-ref
       maxdiff = MaxVal( i )
       If (maxdiff == 0 ) test_result = TEST_PASS
       If (maxdiff /= 0 ) test_result = TEST_FAIL
       Write(test_msg,*) "check array", maxdiff
    End If

    If ( test_result == TEST_PASS )Then
       ! print *, "perm = ", perm
       perm = perm - permref
       maxdiff = MaxVal( perm )
       If (maxdiff == 0 ) test_result = TEST_PASS
       If (maxdiff /= 0 ) test_result = TEST_FAIL
       Write(test_msg,*) "maxval(res-ref) =", maxdiff
    End If

  End Subroutine XMPH_test_sm_qsort

  ! [+] routine : XMPH_test_sm_assemble ---------------------------------------------
  !
  !> test routine sparse_matrix_assemble
  !!
  !! check that spA1, with 3 added entries
  !! (4,3,10.) (4,3,5.) (4,2,7.) is rigthly assembled.
  !! 
  Subroutine XMPH_test_sm_assemble
    
    Integer    :: iinfo
    MPH_INT :: k
    MPH_INT :: nnz
    Real       :: eps = 1.e-15
    XMPH_FLOAT :: ftmp, ftmp2, ftmp3
    Real       :: diff

    Type(XMPH_sparse_matrix_t) :: sm
    Type(XMPH_sparse_matrix_t) :: smref


    ! Setup
    nnz = XMPH_spA1%nnz
    test_result = TEST_SKIP

#if XMPH_HAVE_ARITH_D 
    ftmp  = 10.d0
    ftmp2 = 5.d0
    ftmp3 = 7.d0
#elif XMPH_HAVE_ARITH_S
    ftmp = 10.
    ftmp2 = 5.
    ftmp3 = 7.
#elif XMPH_HAVE_ARITH_Z
    ftmp = CMPLX(10.d0,0.d0 )
    ftmp2 = CMPLX(5.d0,0.d0 )
    ftmp3 = CMPLX(7.d0,0.d0 )
#elif XMPH_HAVE_ARITH_C
    ftmp = CMPLX(10.,0.)
    ftmp2 = CMPLX(5.,0.)
    ftmp3 = CMPLX(7.,0.)
#endif

    ! sm
    Call XMPH_sm_dup(sm,XMPH_spA1,iinfo)
    If (iinfo < 0) Return
    Call XMPH_sm_realloc(sm,nnz+3,iinfo)
    If (iinfo < 0) Return

    Call XMPH_sm_insertEntry(sm,nnz+1, 4, 3, ftmp, iinfo )
    Call XMPH_sm_insertEntry(sm,nnz+2, 4, 3, ftmp2, iinfo )
    Call XMPH_sm_insertEntry(sm,nnz+3, 4, 2, ftmp3, iinfo )
    If (iinfo < 0) Return

    ! smref
    Call XMPH_sm_dup(smref,XMPH_spA1,iinfo)
    If (iinfo < 0) Return
    Do k= 1, nnz
       If ( (smref%i(k) == 4 ).and.( smref%j(k) == 3)) Then 
          smref%v(k) = smref%v(k) + ftmp + ftmp2
       End If
       If ( (smref%i(k) == 4 ).and.( smref%j(k) == 2)) Then 
          smref%v(k) = smref%v(k) + ftmp3
       End If
    End Do
    test_result = TEST_UNSET

    ! Run
    ! Call XMPH_sm_mmwrite(sm,0,iinfo)
    Call XMPH_sm_assemble(sm,sm%n,iinfo)

    ! Check
    Call XMPH_sm_convert(smref,SM_FMT_CSR,iinfo)
    ! Call XMPH_sm_mmwrite(sm,0,iinfo)
    ! Write(0,*) "nnz=", sm%nnz
    ! Write(0,*) "i=", sm%i
    ! Write(0,*) "j=", sm%j
    ! Write(0,*) "v=", sm%v
    ! Write(0,*) "csr=", sm%csr
    ! Call XMPH_sm_mmwrite(smref,0,iinfo)
    ! Write(0,*) "csr=", smref%csr

    test_result = TEST_PASS
    sm%i(1:nnz) = sm%i(1:nnz)-smref%i(1:nnz)
    sm%j(1:nnz) = sm%j(1:nnz)-smref%j(1:nnz)
    sm%v(1:nnz) = sm%v(1:nnz)-smref%v(1:nnz)
    k= smref%m+1
    sm%csr(1:k) = sm%csr(1:k) -smref%csr(1:k)
    diff = XMPH_norm2_diff(sm%v,sm%v,nnz)

    If ( Maxval(sm%i(1:nnz))   /= 0         ) test_result = TEST_FAIL 
    If ( Maxval(sm%i(1:nnz))   /= 0         ) test_msg    = "wrong rows"

    If ( Maxval(sm%j(1:nnz))   /= 0         ) test_result = TEST_FAIL 
    If ( Maxval(sm%j(1:nnz))   /= 0         ) test_msg    = "wrong cols"

    If ( diff > eps  ) test_result = TEST_FAIL 
    If ( diff > eps  ) test_msg    = "wrong vals"

    If ( Maxval(sm%csr(1:k)) /= 0 ) test_result = TEST_FAIL 
    If ( Maxval(sm%csr(1:k)) /= 0 ) test_msg    = "wrong csr"

    If ( nnz   /= smref%nnz       ) test_result = TEST_FAIL 
    If ( nnz   /= smref%nnz       ) test_msg    = "wrong nnz"


  End Subroutine XMPH_test_sm_assemble




  ! [+] routine : XMPH_test_sm_symStruct --------------------------------------------
  !
  !> test routine sparse_matrix_symStruct
  !!
  !! check that spA1, with and added entry (3,4,0.) is rigthly symmetrised.
  !! 
  Subroutine XMPH_test_sm_symStruct
    

    Integer    :: iinfo
    MPH_INT :: nnz
    Real       :: eps = 1.e-15
    Type(XMPH_sparse_matrix_t) :: sm
    Type(XMPH_sparse_matrix_t) :: smref
    Real       :: diff

    ! Setup
    nnz = XMPH_spA1%nnz
    test_result = TEST_SKIP
    
    ! sm
    Call XMPH_sm_dup(sm,XMPH_spA1,iinfo)
    If (iinfo < 0) Return
    Call XMPH_sm_realloc(sm,nnz+1,iinfo)
    If (iinfo < 0) Return
    Call XMPH_sm_insertEntry(sm,nnz+1, 3, 4, XMPH_FLOATZERO, iinfo )
    If (iinfo < 0) Return
    ! Call XMPH_sm_convert(sm,SM_FMT_CSR,iinfo) ! yoh : to test in csr 

    ! smref
    Call XMPH_sm_dup(smref,XMPH_spA1,iinfo)
    If (iinfo < 0) Return
    Call XMPH_sm_realloc(smref,nnz+3,iinfo)
    If (iinfo < 0) Return
    Call XMPH_sm_insertEntry(smref,nnz+1, 3, 4, XMPH_FLOATZERO, iinfo )
    If (iinfo < 0) Return
    Call XMPH_sm_insertEntry(smref,nnz+2, 1, 3, XMPH_FLOATZERO, iinfo )
    If (iinfo < 0) Return
    Call XMPH_sm_insertEntry(smref,nnz+3, 2, 4, XMPH_FLOATZERO, iinfo )
    If (iinfo < 0) Return
    test_result = TEST_UNSET

    ! Run
    Call XMPH_sm_symStruct(sm,iinfo)

    ! Check
    Call XMPH_sm_convert(sm   ,SM_FMT_CSC,iinfo)
    Call XMPH_sm_convert(smref,SM_FMT_CSC,iinfo)
    ! Call XMPH_sm_mmwrite(sm,6,iinfo)
    ! Call XMPH_sm_mmwrite(smref,6,iinfo)

    test_result = TEST_PASS
    sm%i(1:nnz) = sm%i(1:nnz)-smref%i(1:nnz)
    sm%j(1:nnz) = sm%j(1:nnz)-smref%j(1:nnz)
    sm%v(1:nnz) = sm%v(1:nnz)-smref%v(1:nnz)

    diff = XMPH_norm2_diff(sm%v,sm%v,nnz)

    If ( Maxval(sm%i(1:nnz)) /= 0 ) test_msg    = "wrong rows"
    If ( Maxval(sm%j(1:nnz)) /= 0 ) test_result = TEST_FAIL 
    If ( Maxval(sm%j(1:nnz)) /= 0 ) test_msg    = "wrong cols"
    If ( diff > eps               ) test_result = TEST_FAIL 
    If ( diff > eps               ) test_msg    = "wrong vals"

  End Subroutine XMPH_test_sm_symStruct


  ! [+] routine : XMPH_test_sm_matvect_unsym ------------------------------------------
  !
  !> testing the matrix vector product : A.x = y, with A unsymmetric 
  Subroutine XMPH_test_sm_matvect_unsym

    Call XMPH_test_sm_matvect(XMPH_spA1,XMPH_x1,XMPH_b1,1.e-6)

  End Subroutine XMPH_test_sm_matvect_unsym


  ! [+] routine : XMPH_test_sm_matvect_sym -------------------------------------------
  !
  !> testing the matrix vector product : A.x = y, with A symmetric 
  Subroutine XMPH_test_sm_matvect_sym

    Call XMPH_test_sm_matvect(XMPH_spA2,XMPH_x2,XMPH_b2,1.e-6)

  End Subroutine XMPH_test_sm_matvect_sym


  ! [+] routine : XMPH_test_sm_matvect_spd -------------------------------------------
  !
  !> testing the matrix vector product : A.x = y, with A spd
  Subroutine XMPH_test_sm_matvect_spd

    Call XMPH_test_sm_matvect(XMPH_spA3,XMPH_x3,XMPH_b3,1.e-6)

  End Subroutine XMPH_test_sm_matvect_spd


  ! [+] routine : XMPH_test_sm_matvect --------------------------------------------
  !
  !> test that sparce matrix / vector product "A.x = b"
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
  Subroutine XMPH_test_sm_matvect( spA, x, b, eps )

    !* Arguments *!
    Type(XMPH_sparse_matrix_t), Intent(in) :: spA
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

    Call XMPH_sm_vectorproduct(spA,0,0, x, y, status )
    If ( status /= 0 ) Return

    ! Write (*,*)
    ! Write (*,*) "b=", b%v
    ! Write (*,*) "y=", y%v
    ! Write (*,*) "y-b=", y%v-b%v
    ! Write (*,*)
    Write(test_msg,*) " | y-b | =",&
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


  End Subroutine XMPH_test_sm_matvect


End Module XMPH_test_sparse_matrix_mod
