! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"

!> module which shares a set of fixtures for the tests.
!! namely small default linear systems
!!
!!----
!!
!! @par system 1 : unsymmetric linear system
!!      
!! @verbatim
!!       - REAL -         |         - COMPLEX -
!!  1 0 0 0 0     1    1  |   1   0   0   0   0     1+i      1 +i    
!!  0 3 0 0 0     1    3  |   0   3   0   0   0     1+i      3 +3i
!!  2 0 5 0 0  .  1 =  7  |  2+i  0   5   0   0  .  1+i   =  6 +8i
!!  0 4 6 7 0     1   17  |   0  4+i 6+i  7   0     1+i      15+7i
!!  0 0 0 0 8     1    8  |   0   0   0   0   8     1+i      8 +8i
!!
!! Slist=(3,4,5)
!!
!!     5  0  0            |   5  0  0
!! S = 6  7  0            |  6+i 7  0 
!!     0  0  8            |   0  0  8      
!!
!! @endverbatim
!!
!!
!!----
!!
!! @par system 2 : symmetric linear system
!!
!! @verbatim
!!       - REAL -         |         - COMPLEX -
!!  1 0 2 0 0     1    3  |   1   0  2+i  0   0     1+i      2 +4i      
!!  0 3 0 4 0     1    7  |   0   3   0  4+i  0     1+i      6 +8i        
!!  2 0 5 6 0  .  1 = 13  |  2+i  0   5  6+i  0  .  1+i   =  11+15i
!!  0 4 6 7 0     1   17  |   0  4+i 6+i  7   0     1+i      15+19i
!!  0 0 0 0 8     1    8  |   0   0   0   0   8     1+i      8 +8i  
!!
!! Slist=(3,4,5)
!!
!!     1  6  0            |  2-4i     6+i  0
!! S = 6 5/3 0            |   6+i  2-8/3i  0 
!!     0  0  8            |   0    0       8      
!!
!! @endverbatim
!!
!!----
!!
!! @par system 3 : real symmetric positive definite linear system (poisson)
!!
!! @verbatim
!!               --  Linear system  --
!!       - REAL -               |  - COMPLEX is a copy of the REAL case -
!!   2 -1  0  0  0     1    1                
!!  -1  2 -1  0  0     1    0          
!!   0 -1  2 -1  0  .  1 =  0
!!   0  0 -1  2 -1     1    0
!!   0  0  0 -1  2     1    1    
!!
!!               --  Schur complement  --
!!  Slist = (3,4,5)
!!
!!  S = A22 - A21.A11^{-1}.A12                             
!!
!!  A11=A(1:2,1:2) A12=A(1:2,3:5)
!!  A21=A(3:5,1:2) A22=A(3:5,3:5)
!!  A11^{-1} = 1/3 * | 2  1 |
!!                   | 1  2 |
!!       4/3 -1  0  
!!  S =   -1  2 -1
!!        0  -1  2
!!
!! @endverbatim
!!
Module XMPH_test_fixture_mod

  !* Module(s) used *!
  Use XMPH_dense_matrix_mod
  Use XMPH_sparse_matrix_mod

  !* No implicit typing *!
  Implicit none

  !* Shared variables *!

  Type(XMPH_dense_matrix_t), Public, Save :: XMPH_A1  
  Type(XMPH_sparse_matrix_t), Public, Save :: XMPH_spA1  
  Type(XMPH_dense_matrix_t), Public, Save :: XMPH_x1 
  Type(XMPH_dense_matrix_t), Public, Save :: XMPH_b1 
  Integer, Public, Save    :: XMPH_nS1 
  Integer, Public, Pointer :: XMPH_S1list(:)
  Type(XMPH_dense_matrix_t), Public, Save :: XMPH_S1   
  Type(XMPH_sparse_matrix_t), Public, Save :: XMPH_spS1  

  Type(XMPH_dense_matrix_t), Public, Save :: XMPH_A2  
  Type(XMPH_sparse_matrix_t), Public, Save :: XMPH_spA2  
  Type(XMPH_dense_matrix_t), Public, Save :: XMPH_x2 
  Type(XMPH_dense_matrix_t), Public, Save :: XMPH_b2 
  Integer, Public, Save    :: XMPH_nS2 
  Integer, Public, Pointer :: XMPH_S2list(:)
  Type(XMPH_dense_matrix_t), Public, Save :: XMPH_S2   
  Type(XMPH_sparse_matrix_t), Public, Save :: XMPH_spS2  

  Type(XMPH_dense_matrix_t), Public, Save :: XMPH_A3 
  Type(XMPH_sparse_matrix_t), Public, Save :: XMPH_spA3
  Type(XMPH_dense_matrix_t), Public, Save :: XMPH_x3 
  Type(XMPH_dense_matrix_t), Public, Save :: XMPH_b3 
  Integer, Public, Save    :: XMPH_nS3 
  Integer, Public, Pointer :: XMPH_S3list(:)
  Type(XMPH_dense_matrix_t), Public, Save :: XMPH_S3   
  Type(XMPH_sparse_matrix_t), Public, Save :: XMPH_spS3

  !* Private variables *!

  Integer, Private, Save    :: init = 0
  Integer, Private, Pointer :: segfault (:)
  Integer, Private, Parameter :: wp = XMPH_FLOATKIND

  !* Access specifiers to routines *!
  Public :: XMPH_test_fixture_init
  Public :: XMPH_test_fixture_exit
  Public :: XMPH_norm_inf

Contains

  ! [+] routine : XMPH_test_fixture_init ------------------------------------
  !
  !> Initialize the test_fixture module
  !! 
  Subroutine XMPH_test_fixture_init

    Integer :: status
    Integer, Parameter :: n  = 5 ! order of the linear systems
    Integer, Parameter :: nz = 8 ! number of entries in the matrices A1 & A2

    ! End of header ------------------------------------------------------------

    ! exit early if module is already initialized
    If ( init > 0 )Then
       init = init+1
       Return
    End If

    !-------------------------------------------------------------------------
    ! [1] Set the unsymmetric linear system
    !-------------------------------------------------------------------------

    ! 
    Call XMPH_dm_create( XMPH_x1, n, 1, n, status )
    If ( status /= 0 ) Goto 9999
    Call XMPH_dm_create( XMPH_b1, n, 1, n, status )
    If ( status /= 0 ) Goto 9999

    Call XMPH_dm_create( XMPH_A1, n, n, n, status )
    If ( status /= 0 ) Goto 9999
    Call XMPH_sm_create( XMPH_spA1, nz, status )
    If ( status /= 0 ) Goto 9999

#if (XMPH_HAVE_ARITH_D) || (XMPH_HAVE_ARITH_S)
    ! REAL
    XMPH_x1%v = (/ 1 , 1 , 1 , 1 , 1 /)
    XMPH_b1%v = (/ 1 , 3 , 7 , 17 , 8 /)
    XMPH_A1%v = (/     &
         1, 0, 2, 0, 0 ,&  ! col 1
         0, 3, 0, 4, 0 ,&  ! col 2
         0, 0, 5, 6, 0 ,&  ! col 3
         0, 0, 0, 7, 0 ,&  ! col 4
         0, 0, 0, 0, 8 /)  ! col 5

    XMPH_spA1%m     = n
    XMPH_spA1%n     = n
    XMPH_spA1%i(:)  =(/1, 3, 2, 4, 3, 4, 4, 5/)
    XMPH_spA1%j(:)  =(/1, 1, 2, 2, 3, 3, 4, 5/)
    XMPH_spA1%v(:)  =(/1, 2, 3, 4, 5, 6, 7, 8/)

    XMPH_nS1 = 3
    Allocate(XMPH_S1list(XMPH_nS1))
    XMPH_S1list = (/ 3, 4 , 5 /)
    Call XMPH_dm_create(XMPH_S1,XMPH_nS1,XMPH_nS1,XMPH_nS1,status)
    XMPH_S1%v = (/ &
         5, 6, 0 ,&  ! col 1
         0, 7, 0 ,&  ! col 2
         0, 0, 8 /)  ! col 3

    Call XMPH_sm_create(XMPH_spS1,4,status)
    XMPH_SpS1%n = XMPH_nS1
    XMPH_SpS1%i = (/ 1, 2, 3, 2 /) 
    XMPH_SpS1%j = (/ 1, 2, 3, 1 /)
    XMPH_SpS1%v = (/ 5, 7, 8, 6 /)

#elif (XMPH_HAVE_ARITH_C) || (XMPH_HAVE_ARITH_Z)

    ! COMPLEX
    XMPH_x1%v = (/ (1,1),(1,1),(1,1),(1,1),(1,1)  /)
    XMPH_b1%v = (/ (1,1),(3,3),(6,8),(15,19),(8,8) /)
    XMPH_A1%v = (/     &
         (1,0), (0,0), (2,1), (0,0), (0,0),&  ! col 1
         (0,0), (3,0), (0,0), (4,1), (0,0),&  ! col 2
         (0,0), (0,0), (5,0), (6,1), (0,0),&  ! col 3
         (0,0), (0,0), (0,0), (7,0), (0,0),&  ! col 4
         (0,0), (0,0), (0,0), (0,0), (8,0)/)  ! col 5

    XMPH_spA1%m     = n
    XMPH_spA1%n     = n
    XMPH_spA1%i(:)  =(/1, 3, 2, 4, 3, 4, 4, 5/)
    XMPH_spA1%j(:)  =(/1, 1, 2, 2, 3, 3, 4, 5/)
    XMPH_spA1%v(:)  =(/(1,0), (2,1), (3,0), (4,1), (5,0), (6,1), (7,0), (8,0)/)

    XMPH_nS1 = 3
    Allocate(XMPH_S1list(XMPH_nS1))
    XMPH_S1list = (/ 3, 4 , 5 /)
    Call XMPH_dm_create(XMPH_S1,XMPH_nS1,XMPH_nS1,XMPH_nS1,status)
    XMPH_S1%v = (/ &
         (5,0), (6,1), (0,0),&  ! col 1
         (0,0), (7,0), (0,0),&  ! col 2
         (0,0), (0,0), (8,0)/)  ! col 3

    Call XMPH_sm_create(XMPH_spS1,4,status)
    XMPH_SpS1%n = XMPH_nS1
    XMPH_SpS1%i = (/ 1, 2, 3, 2 /) 
    XMPH_SpS1%j = (/ 1, 2, 3, 1 /)
    XMPH_SpS1%v = (/ (5,0), (7,0), (8,0), (6,1) /)

#endif

    !-------------------------------------------------------------------------
    ! [2] symmetric system
    !-------------------------------------------------------------------------

    ! 
    Call XMPH_dm_create( XMPH_x2, n, 1, n, status )
    If ( status /= 0 ) Goto 9999
    Call XMPH_dm_create( XMPH_b2, n, 1, n, status )
    If ( status /= 0 ) Goto 9999

    Call XMPH_dm_create( XMPH_A2, n, n, n, status )
    If ( status /= 0 ) Goto 9999
    Call XMPH_sm_create( XMPH_spA2, nz , status )
    If ( status /= 0 ) Goto 9999


#if ( XMPH_HAVE_ARITH_D ) || ( XMPH_HAVE_ARITH_S )
    ! REAL

    !! @verbatim
    !!  1 0 2 0 0     1    3       
    !!  0 3 0 4 0     1    7          
    !!  2 0 5 6 0  .  1 = 13
    !!  0 4 6 7 0     1   17
    !!  0 0 0 0 8     1    8    
    !! @endverbatim

    XMPH_x2%v = (/ 1 , 1 , 1 , 1 , 1 /)
    XMPH_b2%v = (/ 3 , 7 ,13 ,17 , 8 /)

    XMPH_A2%sym    = DM_SYM_IsSymmetric
    XMPH_A2%stored = DM_STORED_LOWER
    XMPH_A2%v = (/     &
         1, 0, 2, 0, 0 ,&  ! col 1
         0, 3, 0, 4, 0 ,&  ! col 2
         0, 0, 5, 6, 0 ,&  ! col 3
         0, 0, 0, 7, 0 ,&  ! col 4
         0, 0, 0, 0, 8 /)  ! col 5

    XMPH_spA2%sym   = SM_SYM_IsSymmetric
    XMPH_spA2%m     = n
    XMPH_spA2%n     = n
    XMPH_spA2%i(:)  =(/1, 3, 2, 4, 3, 4, 4, 5/)
    XMPH_spA2%j(:)  =(/1, 1, 2, 2, 3, 3, 4, 5/)
    XMPH_spA2%v(:)  =(/1, 2, 3, 4, 5, 6, 7, 8/)

    XMPH_nS2 = 3
    Allocate(XMPH_S2list(XMPH_nS2))
    XMPH_S2list = (/ 3, 4 , 5 /)
    Call XMPH_dm_create(XMPH_S2,XMPH_nS2,XMPH_nS2,XMPH_nS2,status)
    XMPH_S2%v = (/ &
         1._wp,        6._wp  ,  0._wp ,& ! col 1
         6._wp,    5._wp/3._wp,  0._wp ,& ! col 2
         0._wp,        0._wp  ,  8._wp /) ! col 3

#elif (XMPH_HAVE_ARITH_C) || (XMPH_HAVE_ARITH_Z)

    ! COMPLEX
    XMPH_x2%v = (/ (1,1),(1,1),(1,1),(1,1),(1,1)  /)
    XMPH_b2%v = (/ (2,4),(6,8),(11,15),(15,19),(8,8) /)

    XMPH_A2%sym    = DM_SYM_IsSymmetric
    XMPH_A2%stored = DM_STORED_LOWER
    XMPH_A2%v = (/     &
         (1,0), (0,0), (2,1), (0,0), (0,0),&  ! col 1
         (0,0), (3,0), (0,0), (4,1), (0,0),&  ! col 2
         (0,0), (0,0), (5,0), (6,1), (0,0),&  ! col 3
         (0,0), (0,0), (0,0), (7,0), (0,0),&  ! col 4
         (0,0), (0,0), (0,0), (0,0), (8,0)/)  ! col 5

    XMPH_spA2%sym   = SM_SYM_IsSymmetric
    XMPH_spA2%m     = n
    XMPH_spA2%n     = n
    XMPH_spA2%i(:)  =(/1, 3, 2, 4, 3, 4, 4, 5/)
    XMPH_spA2%j(:)  =(/1, 1, 2, 2, 3, 3, 4, 5/)
    XMPH_spA2%v(:)  =(/(1,0), (2,1), (3,0), (4,1), (5,0), (6,1), (7,0), (8,0)/)

    XMPH_nS2 = 3
    Allocate(XMPH_S2list(XMPH_nS2))
    XMPH_S2list = (/ 3, 4 , 5 /)
    Call XMPH_dm_create(XMPH_S2,XMPH_nS2,XMPH_nS2,XMPH_nS2,status)
    ! col 1
    XMPH_S2%v(1) = CMPLX(+2._wp,-4._wp,KIND=wp)
    XMPH_S2%v(2) = CMPLX(+6._wp,+1._wp,KIND=wp)
    XMPH_S2%v(3) = CMPLX(+0._wp,+0._wp,KIND=wp) 
    ! col 2
    XMPH_S2%v(4) = CMPLX(+6._wp,+1._wp,KIND=wp) 
    XMPH_S2%v(5) = CMPLX(+2._wp,-8._wp/3._wp,KIND=wp)    
    XMPH_S2%v(6) = CMPLX(+0._wp,+0._wp,KIND=wp)
    ! col 3
    XMPH_S2%v(7) = CMPLX(+0._wp,+0._wp,KIND=wp) 
    XMPH_S2%v(8) = CMPLX(+0._wp,+0._wp,KIND=wp) 
    XMPH_S2%v(9) = CMPLX(+8._wp,+0._wp,KIND=wp)

#endif

    !-------------------------------------------------------------------------
    ! [3] SPD system
    !-------------------------------------------------------------------------

    ! 
    Call XMPH_dm_create( XMPH_x3, n, 1, n, status )
    If ( status /= 0 ) Goto 9999
    Call XMPH_dm_create( XMPH_b3, n, 1, n, status )
    If ( status /= 0 ) Goto 9999

    Call XMPH_dm_create( XMPH_A3, n, n, n, status )
    If ( status /= 0 ) Goto 9999
    Call XMPH_sm_create( XMPH_spA3, 9 , status )
    If ( status /= 0 ) Goto 9999

    !! @verbatim
    !!   2 -1  0  0  0     1    1                
    !!  -1  2 -1  0  0     1    0          
    !!   0 -1  2 -1  0  .  1 =  0
    !!   0  0 -1  2 -1     1    0
    !!   0  0  0 -1  2     1    1    
    !! @endverbatim
    !!

    XMPH_x3%v = (/ 1 , 1 , 1 , 1 , 1 /)
    XMPH_b3%v = (/ 1 , 0 , 0 , 0 , 1 /)

    XMPH_A3%sym    = DM_SYM_IsSPD
    XMPH_A3%stored = DM_STORED_LOWER
    XMPH_A3%v = (/     &
         2,-1, 0, 0, 0 ,&  ! col 1
         0, 2,-1, 0, 0 ,&  ! col 2
         0, 0, 2,-1, 0 ,&  ! col 3
         0, 0, 0, 2,-1 ,&  ! col 4
         0, 0, 0, 0, 2 /)  ! col 5

    XMPH_spA3%sym   = SM_SYM_IsSPD
    XMPH_spA3%m     = n
    XMPH_spA3%n     = n
    XMPH_spA3%i(:)  =(/1, 2, 3, 4, 5,  2, 3, 4, 5/)
    XMPH_spA3%j(:)  =(/1, 2, 3, 4, 5,  1, 2, 3, 4/)
    XMPH_spA3%v(:)  =(/2, 2, 2, 2, 2, -1,-1,-1,-1/)

    XMPH_nS3 = 3
    Allocate(XMPH_S3list(XMPH_nS3))
    XMPH_S3list = (/ 3, 4 , 5 /)
    Call XMPH_dm_create(XMPH_S3,XMPH_nS3,XMPH_nS3,XMPH_nS3,status)
    XMPH_S3%v = (/ & 
         4._wp/3._wp , -1._wp,  0._wp ,& ! col 1
           -1._wp    ,  2._wp, -1._wp ,& ! col 2
            0._wp    , -1._wp,  2._wp /) ! col 3

    !-------------------------------------------------------------------------
    ! [-] Finish
    !-------------------------------------------------------------------------

9999 Continue
    If ( status /= 0 ) segfault( 1 ) = 1

  End Subroutine XMPH_test_fixture_init

  ! [+] routine : test_fixture_exit ------------------------------------
  !
  !> cleanup the test_fixture module
  !! 
  Subroutine XMPH_test_fixture_exit
    
    Integer :: status

    ! decrement module usage
    If ( init > 0 ) init = init - 1

    ! free memory 
    If ( init <= 0 ) Then

       Call XMPH_dm_free( XMPH_x1, status )
       Call XMPH_dm_free( XMPH_x1, status )
       Call XMPH_dm_free( XMPH_A1, status )
       Call XMPH_sm_free( XMPH_spA1, status )

       Call XMPH_dm_free( XMPH_x2, status )
       Call XMPH_dm_free( XMPH_x2, status )
       Call XMPH_dm_free( XMPH_A2, status )
       Call XMPH_sm_free( XMPH_spA2, status )

       Call XMPH_dm_free( XMPH_x3, status )
       Call XMPH_dm_free( XMPH_x3, status )
       Call XMPH_dm_free( XMPH_A3, status )
       Call XMPH_sm_free( XMPH_spA3, status )

    End If

  End Subroutine XMPH_test_fixture_exit


  ! [+] function : XMPH_norm_inf --------------------------------------------------
  !
  !> compute norm inf of 2 arrays
  !!
  !! @param [in] a1   the first  array
  !! @param [in] a1   the second array
  !! @return          the infinite norm of their difference
  !!
  Function XMPH_norm_inf(  a1 , a2, len )

    Real(kind=8) :: XMPH_norm_inf
    XMPH_FLOAT :: a1 (:), a2(:)
    Integer      :: len

    Integer :: i
    Real(kind=8) :: diff
    !
    XMPH_norm_inf = 0.d0

    Do i = 1, len
       diff = Abs( a1(i) - a2(i) )
       If (diff > XMPH_norm_inf) &
            XMPH_norm_inf = diff
    End Do

  End Function XMPH_norm_inf


  ! [+] function : XMPH_norm2_diff --------------------------------------------------
  !
  !> compute norm 2 of the difference of 2 arrays
  !!
  !! @param [in] a1   the first  array
  !! @param [in] a1   the second array
  !! @return          the norm 2 of their difference
  !!
  Function XMPH_norm2_diff(  a1 , a2, len )

    Real(kind=XMPH_FLOATKIND) :: XMPH_norm2_diff
    XMPH_FLOAT , Intent(in) :: a1 (:), a2(:)
    Integer      , Intent(in) :: len

    Integer :: i
    XMPH_FLOAT :: da (len)
    XMPH_FLOAT, External :: XMPH_DOT

    !
    Do i= 1,len
       da(i) = a1(i) - a2(i)
    End Do

    XMPH_norm2_diff = Real(XMPH_DOT(len,da,1,da,1), XMPH_FLOATKIND)
    XMPH_norm2_diff = SQRT(XMPH_norm2_diff)

  End Function XMPH_norm2_diff



End Module XMPH_test_fixture_mod


