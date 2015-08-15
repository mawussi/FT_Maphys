! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"

! [+] module : dds_lapack_mod -----------------------------------------------
!
!> Public interface for a generic dense direct solver
!!
!! @warning Convention on symmetric/SPD matrices with ALL entries.
!! If a symmetric/SPD matrices references ALL the entries,
!! by convention with only use the UPPER entries (uplo='U').
!! 
Module XMPH_dds_lapack_mod

  !* Module used *!
  Use MPH_dds_lapack_enum
  Use XMPH_dense_matrix_mod, Only : &
       XMPH_dense_matrix_t        ! derived type

  !* No Implicit Typing *!

  Implicit none

  !-----------------------------------------------------------------------------
  ! [+++] Defined constants and derived types
  !-----------------------------------------------------------------------------


  !> filename
  Character(len=MAPHYS_STRL), Parameter, Private :: &
       FLNAME= "XMPH_ARITHmph_dds_lapack_mod.F90"

  !-------------------------------------------------------------------
  ! [++] Defined derived types
  !-------------------------------------------------------------------

  ! [+] type : maphys_lapack_solver_t --------------------------------
  !
  Type XMPH_lapack_t; sequence
     
     !> options
     Integer :: option(MAPHYS_LAPACK_OPTION_SIZE)

     !> state (attributes)
     Integer :: state(MAPHYS_LAPACK_STATE_SIZE)

     !> pivots (produced during factorization)
     Integer, Pointer :: IPIV (:)

     !> statistics
     Real(kind=8) :: rinfo(MAPHYS_LAPACK_RINFO_SIZE)

  End Type XMPH_lapack_t

  !-----------------------------------------------------------------------------
  ! [+++] Methods
  !-----------------------------------------------------------------------------

  Public :: XMPH_LAPACK_Init
  Public :: XMPH_LAPACK_Exit
  Public :: XMPH_LAPACK_Get_pivots
  Public :: XMPH_LAPACK_Factorize
  Public :: XMPH_LAPACK_Solve

Contains

  ! [+] routine :  XMPH_LAPACK_Init ------------------------------------------
  !
  !> initialize the instance
  !!
  !! Initialize the instance.
  !!
  !! @param[in,out] dds_lapack   the lapack instance
  !! @param[in,out] info         the status of the routine
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_LAPACK_Init ( & ! intents
       dds_lapack ,               & ! inout
       info                       & ! out
       )

    !* Arguments *!

    Type(XMPH_lapack_t), intent(inout) :: dds_lapack
    integer              , intent(  out) :: info

    !- End of header -------------------------------------------------

    !-----------------------------------------------------------------
    ! [1] Nullify pointers
    !-----------------------------------------------------------------

    Nullify(dds_lapack%ipiv)

    !-----------------------------------------------------------------
    ! [2] Give default values to attributes
    !-----------------------------------------------------------------

    dds_lapack%state(LAPACK_STATE_isValid      ) = MPH_TRUE
    dds_lapack%state(LAPACK_STATE_Factorized   ) = MPH_FALSE

    dds_lapack%option(LAPACK_OPTION_CompStats) = MPH_TRUE

    !-----------------------------------------------------------------
    ! [3] Exit routine
    !-----------------------------------------------------------------

    info = 0
    Return

  End Subroutine XMPH_LAPACK_Init


  ! [+] routine :  XMPH_LAPACK_Exit ------------------------------------------
  !
  !> Delete the instance
  !!
  !! Delete the instance. Simply deallocate the temporary memory used.
  !!
  !! @param[in,out] dds_lapack   the lapack instance
  !! @param[in,out] info         the status of the routine
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_LAPACK_Exit ( & ! intents
       dds_lapack ,               & ! inout
       info                       & ! out
       )
    Implicit None

    !* Subroutine arguments *!
    Type(XMPH_lapack_t), intent(inout) :: dds_lapack
    integer              , intent(  out) :: info

    !- End of header -------------------------------------------------

    !-----------------------------------------------------------------
    ! [1] Deallocate or nullify everything
    !-----------------------------------------------------------------

    If (Associated(dds_lapack%ipiv)) Deallocate(dds_lapack%ipiv)

    !-----------------------------------------------------------------
    ! [2] Exit routine
    !-----------------------------------------------------------------

    info = 0
    Return

  End Subroutine XMPH_LAPACK_Exit

  ! [+] routine :  XMPH_LAPACK_get_pivots ------------------------------------
  !
  !> Get the pivots computed during the factorization.
  !! 
  !! @param[in,out] dds_lapack  the lapack instance
  !! @param[in,out] piv         the pivots
  !! @param[   out] info        the routine status
  !!
  !! @author Yohan Lee-tin-yien
  !!
  subroutine XMPH_LAPACK_get_pivots (dds_lapack, piv, info)

    !* Module *!
    Use mph_error_mod
    Implicit None

    !* Arguments *!

    Type(XMPH_lapack_t), intent(inout) :: dds_lapack
    integer, pointer           , intent(inout) :: piv (:)
    integer                    , intent(  out) :: info
    
    !- End of header -------------------------------------------------

    !-----------------------------------------------------------------
    ! [1] Verify 
    !-----------------------------------------------------------------

    CHCKASSRT( Associated(dds_lapack%ipiv), info )
    If( info < 0) Return

    !-----------------------------------------------------------------
    ! [2] Link with the pivots
    !-----------------------------------------------------------------

    Nullify(piv)
    piv => dds_lapack%iPIV

  End Subroutine XMPH_LAPACK_get_pivots

  ! [+] routine :  XMPH_LAPACK_factorize -------------------------------------
  !> factorize a dense matrix
  !!
  !! @copydoc DDS_Factorize()
  !!
  !!----- 
  !!
  !! @author Yohan Lee-tin-yien
  !!
  !! @par History
  !!
  !! @verbatim
  !!   Date     :Version: Comments
  !! - 20/01/11 : v0.1b : add symmetric matrice support. (Yohan Lee-tin-yien)
  !! - 2010     : v0.1a : first implementation.           (Yohan Lee-tin-yien)
  !!                      It only handle unsymmetric matrices.
  !! @endverbatim
  subroutine XMPH_LAPACK_factorize(dds_lapack, MAT, info)

    !* Module(s) *!
    Use XMPH_dense_matrix_type
    Use mph_error_mod
    Implicit None

    !* Arguments *!

    Type(XMPH_lapack_t), Intent(inout) :: dds_lapack
    Type(XMPH_dense_matrix_t) , Intent(inout) :: MAT
    Integer              , Intent(  out) :: info

    !* Local Variables *!

    ! Scalars
    Logical :: tst
    Integer :: iinfo
    Integer :: istep
    Integer :: MatrixType
    Integer :: m,n, ldA
    Integer :: min_mn

    MPH_INT   :: lwork         ! sizes of workspace "work"
    MPH_INT   :: liwork        ! sizes of workspace "iwork"
    MPH_INT   :: lrwork        ! sizes of workspace "rwork"
    XMPH_FLOAT :: lworkf        ! lwork but saved into a float

    Logical :: ComputeStats     ! compute the statistics or not
    Real(kind=XMPH_FLOATKIND) :: norm1A ! norm1 of matrix A
    Real(kind=XMPH_FLOATKIND) :: rcond1 ! reciprocical condition number

    ! Arrays
    Integer     , Pointer :: iPIV(:)
    XMPH_FLOAT, Pointer :: A   (:)
    XMPH_FLOAT, Pointer :: work  (:)
    MPH_INT  , Pointer :: iwork (:)
    Real(kind=XMPH_FLOATKIND), Pointer :: rwork (:)

    ! String 
    Character*1 :: uplo

    ! Lapack functions
    Real(kind=XMPH_FLOATKIND), External :: XMPH_ARITHlange
    Real(kind=XMPH_FLOATKIND), External :: XMPH_ARITHlansy

    !- End of header -------------------------------------------------

    !-----------------------------------------------------------------
    ! [0] Verify Object state
    !-----------------------------------------------------------------

    tst = (dds_lapack%state(LAPACK_STATE_isValid) == MPH_TRUE )
    CHCKASSRT(tst, info )
    If( info < 0) Goto 9999

    !-----------------------------------------------------------------
    ! [1] Initialize local variables
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! [1.1] Initialize local scalars / Arrays
    !-----------------------------------------------------------------

    ! Nullify pointers
    Nullify(iPIV)
    Nullify(A)
    Nullify(work)
    Nullify(iwork)
    Nullify(rwork)

    ! Attributes
    MatrixType    = MAT%sym
    ComputeStats  = (dds_lapack%option(LAPACK_OPTION_CompStats) == MPH_TRUE)

    ! Data 
    m   =  MAT%m
    n   =  MAT%n
    ldA =  MAT%ld
    A   => MAT%v

    ! Set uplo
    Select Case (MatrixType)
    Case( DM_SYM_IsSPD, DM_SYM_IsSymmetric )

       If ( MAT%stored == DM_STORED_UPPER ) uplo='U'
       If ( MAT%stored == DM_STORED_LOWER ) uplo='L'
       ! convention, see the warning at top of file for details
       If ( MAT%stored == DM_STORED_ALL   ) uplo='U'

    Case Default

       uplo=' ' ! uplo is not needed 

    End Select

    ! Set default workspaces sizes 
    ! negative values means do not allocate their corresponding workspaces 
    lwork = -1
    liwork = -1
    lrwork = -1

    !-----------------------------------------------------------------
    ! [1.3] Allocate memory, workspace, etc.
    !-----------------------------------------------------------------
    
    !-----------------------------------------------------------------
    ! [1.3.1] allocate ipiv
    !-----------------------------------------------------------------

    min_mn = min(m,n)
    Allocate( iPIV( min_mn ) , STAT = info )
    CHCKASSRT( info == 0, info )
    If( info < 0) Goto 9999

    !-----------------------------------------------------------------
    ! [1.3.2] Query workspaces sizes : lwork, liwork, lrwork
    !-----------------------------------------------------------------
    istep = 132

    ! For the norm : {C,D,Z,C}LAN{GE,SY}
    IF ( ComputeStats ) Then
       Select Case (MatrixType)
       Case( DM_SYM_IsGeneral   ); lrwork  = Max( lrwork, m )
       Case( DM_SYM_IsSPD       ); lrwork  = Max( lrwork, n )
       Case( DM_SYM_IsSymmetric ); lrwork  = Max( lrwork, n )
       Case Default              ; Continue 
       End Select
    End IF

    ! For the factorization : {C,D,Z,C}{GE,SY,PO}TRF
    Select Case (MatrixType)
    Case( DM_SYM_IsGeneral   )
       Continue 
    Case( DM_SYM_IsSPD       )
       Continue
    Case( DM_SYM_IsSymmetric )
       lwork = -1 ! ask for optimal lwork
       Call XMPH_ARITHsytrf( uplo, n, A, ldA, iPIV, lworkf , LWORK, info )
       CHCKASSRT( info == 0, info )
       If( info < 0) Goto 9999
       lwork = Int(Real(lworkf))
    Case Default
       Continue
    End Select

    ! For the condition number : {C,D,Z,C}{GE,SY,PO}CON
    IF ( ComputeStats ) Then

#if (XMPH_HAVE_ARITH_D) || (XMPH_HAVE_ARITH_S)
       Select Case (MatrixType)
       Case( DM_SYM_IsGeneral ) ! {C,D,Z,C}GECON
          lwork  = Max( lwork , 4*n )
          liwork = Max( liwork, n   )
       Case( DM_SYM_IsSPD       ) ! {C,D,Z,C}POCON
          lwork  = Max( lwork , 3*n )
          liwork = Max( liwork, n   )
       Case( DM_SYM_IsSymmetric ) ! {C,D,Z,C}SYCON
          lwork  = Max( lwork , 2*n )
          liwork = Max( liwork, n   )
       Case Default
          Continue 
       End Select
#elif (XMPH_HAVE_ARITH_C) || (XMPH_HAVE_ARITH_Z)
       Select Case (MatrixType)
       Case( DM_SYM_IsGeneral ) ! {C,D,Z,C}GECON
          lwork  = Max( lwork , 2*n )
          lrwork = Max( lrwork, 2*n )
       Case( DM_SYM_IsSPD       ) ! {C,D,Z,C}POCON
          lwork  = Max( lwork , 2*n )
          lrwork = Max( lrwork, n   )
       Case( DM_SYM_IsSymmetric ) ! {C,D,Z,C}SYCON
          lwork  = Max( lwork , 2*n )
       Case Default
          Continue 
       End Select
#else
#error "bad arithmetic"
#endif

    End IF


    !-----------------------------------------------------------------
    ! [1.3.3] Allocate necessary workspaces
    !-----------------------------------------------------------------

    istep = 133
    If ( lwork  > 0 ) Allocate(  work(lwork ), STAT= info )
    CHCKASSRT( info == 0, info )
    If( info < 0) Goto 9999

    If ( liwork > 0 ) Allocate( iwork(liwork), STAT= info )
    CHCKASSRT( info == 0, info )
    If( info < 0) Goto 9999

    If ( lrwork > 0 ) Allocate( rwork(lrwork), STAT= info )
    CHCKASSRT( info == 0, info )
    If( info < 0) Goto 9999

    !-----------------------------------------------------------------
    ! [1.4] Compute statistics : the norm of input matrix
    !-----------------------------------------------------------------

    !> @warning why "rwork" and not "work" here ?
    !! In Lapack documentation, the workspace of {C,D,Z,C}LAN{GE,SY}
    !! is a real in single or double precision workspace.
    If ( ComputeStats ) Then
       istep = 14

       Select Case (MatrixType)
       Case( DM_SYM_IsGeneral )
          norm1A = XMPH_ARITHlange( '1', m, n, A , ldA , rwork )
       Case( DM_SYM_IsSPD       )
          norm1A = XMPH_ARITHlansy('1', uplo, N, A, ldA, rwork )
       Case( DM_SYM_IsSymmetric )
          norm1A = XMPH_ARITHlansy('1', uplo, N, A, ldA, rwork )
       Case Default
          CHCKASSRT( .False. , info )
          If( info < 0) Goto 9999
       End Select

    End If
    
    !-----------------------------------------------------------------
    ! [2] Call the right factorization routine
    !-----------------------------------------------------------------
    istep = 2

    iinfo = 0
    Select Case (MatrixType)
    Case( DM_SYM_IsGeneral )

       Call XMPH_ARITHgetrf ( m, n, A, ldA, iPIV , iinfo )

    Case( DM_SYM_IsSPD       )

       Call XMPH_ARITHpotrf ( uplo, n, A, ldA, iinfo )

    Case( DM_SYM_IsSymmetric )

       Call XMPH_ARITHSYTRF( uplo, n, A, ldA, iPIV, work, lwork, iinfo )

    Case Default

       iinfo = -1

    End Select

    CHCKASSRT( iinfo >= 0, info )
    If( info < 0) Goto 9999

    !-----------------------------------------------------------------
    ! [3] Compute statistics : the reciprocical of the condition number
    !-----------------------------------------------------------------
    istep= 3

    IF ( ComputeStats ) Then
       Select Case (MatrixType)

       Case( DM_SYM_IsGeneral )! call different Fgecon

          !-----------------------------------------------------------------
          ! [3.1] General matrices 
          !-----------------------------------------------------------------
          istep= 31
#if XMPH_HAVE_ARITH_D
          Call Dgecon( '1', n, A, ldA , norm1A, RCOND1, work, iwork, iinfo)
#elif XMPH_HAVE_ARITH_S
          Call Sgecon( '1', n, A, ldA , norm1A, RCOND1, work, iwork, iinfo)
#elif XMPH_HAVE_ARITH_C
          Call Cgecon( '1', n, A, ldA , norm1A, RCOND1, work, rwork, iinfo)
#elif XMPH_HAVE_ARITH_Z
          Call Zgecon( '1', n, A, ldA , norm1A, RCOND1, work, rwork, iinfo)
#else
          iinfo = -istep
#endif
          CHCKASSRT( iinfo >= 0, info )
          If( info < 0) Goto 9999

       Case( DM_SYM_IsSPD       )

          !-----------------------------------------------------------------
          ! [3.2] SPD matrices 
          !-----------------------------------------------------------------

          istep = 32

#if XMPH_HAVE_ARITH_D
          Call Dpocon( uplo, n, A, ldA, norm1a, RCOND1, work, iwork, iinfo)
#elif XMPH_HAVE_ARITH_S
          Call Spocon( uplo, n, A, ldA, norm1a, RCOND1, work, iwork, iinfo)
#elif XMPH_HAVE_ARITH_C
          Call Cpocon( uplo, n, A, ldA, norm1a, RCOND1, work, rwork, iinfo)
#elif XMPH_HAVE_ARITH_Z
          Call Zpocon( uplo, N, A, ldA, norm1a, RCOND1, work, rwork, iinfo)
#else
          iinfo = -istep
#endif
          If( iinfo < 0) Goto 9999 

       Case( DM_SYM_IsSymmetric) ! call different dsycon

          !-----------------------------------------------------------------
          ! [3.3] Symmetric matrices 
          !-----------------------------------------------------------------

          istep = 33
#if XMPH_HAVE_ARITH_D
          Call Dsycon( uplo, n, A, ldA, ipiv, norm1a, RCOND1, work, iwork, iinfo)
#elif XMPH_HAVE_ARITH_S
          Call Ssycon( uplo, n, A, ldA, ipiv, norm1a, RCOND1, work, iwork, iinfo)
#elif XMPH_HAVE_ARITH_C
          Call Csycon( uplo, n, A, ldA, ipiv, norm1a, RCOND1, work, iinfo)
#elif XMPH_HAVE_ARITH_Z
          Call Zsycon( uplo, N, A, ldA, ipiv, norm1a, RCOND1, work, iinfo)
#else
          iinfo = -istep
#endif

          CHCKASSRT( iinfo >= 0, info )
          If( info < 0) Goto 9999

       Case Default

          !-----------------------------------------------------------------
          ! [3.4] Other matrices 
          !-----------------------------------------------------------------

          CHCKASSRT(.False., info )
          If( info < 0) Goto 9999

       End Select

    Endif


    !-----------------------------------------------------------------
    ! [4] Exit routine 
    !-----------------------------------------------------------------

    ! Save data
    dds_lapack%ipiv => iPIV
    IF ( ComputeStats ) Then
       dds_lapack%rinfo(LAPACK_RINFO_Norm1 ) = norm1A
       dds_lapack%rinfo(LAPACK_RINFO_RCond1) = rcond1
    End IF

    ! Update status
    info = 0
    dds_lapack%state(LAPACK_STATE_Factorized ) = MPH_TRUE 

9999 Continue
    ! Free memory
    If(Associated(iwork)) Deallocate(iwork)
    If(Associated(rwork)) Deallocate(rwork)
    If(Associated(work)) Deallocate(work)

  End Subroutine XMPH_LAPACK_factorize


  ! [+] routine :  XMPH_LAPACK_solve -----------------------------------------
  !> solve the dense linear system
  !!
  !! Solve the dense linear system MAT. SOL = RHS.
  !! with Mat previously factorized by XMPH_LAPACK_Factorize().
  !!
  !!-----
  !! 
  !! @param[in,out] dds_lapack 
  !!
  !!   the lapack instance, which contains several data such as the pivots.
  !!
  !! @param[in,out] MAT
  !!
  !!   the matrix previously factorized by XMPH_LAPACK_Factorize()
  !!
  !! @param[in,out] RHS
  !!
  !!   - On input, holds the right-hand-side
  !!   - On output, holds the solution 
  !!
  !! @param[   out] info
  !!
  !!   The routine status
  !!
  !!-----
  !! 
  !! @author Yohan Lee-tin-yien
  !!
  !! @par History
  !!
  !! @verbatim
  !!   Date     :Version: Comments
  !! - 20/01/11 : v0.1b : add symmetric matrices support. (Yohan Lee-tin-yien)
  !! - 2010     : v0.1a : first implementation.           (Yohan Lee-tin-yien)
  !!                      It only handles unsymmetric matrices.
  !! @endverbatim
  subroutine XMPH_LAPACK_solve(dds_lapack, MAT, RHS, info)
    
    !* Modules *!
    Use XMPH_dense_matrix_mod, Only : &
         XMPH_dense_matrix_t, & ! structures
         DM_SYM_IsGeneral, & ! enumerations
         DM_SYM_IsSPD, & 
         DM_SYM_IsSymmetric, &
         DM_STORED_ALL  , &
         DM_STORED_UPPER, &
         DM_STORED_LOWER
    Use mph_error_mod
    Implicit None

    !* Subroutine Arguments *!

    Type(XMPH_lapack_t) , intent(inout) :: dds_lapack
    Type(XMPH_dense_matrix_t)  , Intent(inout) :: MAT
    Type(XMPH_dense_matrix_t)  , Intent(inout) :: RHS
    integer               , intent(  out) :: info

    !* Local Variables *!

    ! Scalars
    Logical :: tst
    Integer :: sym
    Integer :: n, ldA
    Integer :: nrhs, ldb

    ! Arrays
    XMPH_FLOAT, Pointer :: A_1  (:)
    Integer     , Pointer :: ipiv (:)
    XMPH_FLOAT, Pointer :: b    (:)

    ! String 
    Character*1 :: TRANS
    Character*1 :: UPLO

    !- End of header -------------------------------------------------

    !-----------------------------------------------------------------
    ! [0] Verify Object state
    !-----------------------------------------------------------------

    tst = dds_lapack%state(LAPACK_STATE_isValid ) == MPH_TRUE
    CHCKASSRT( tst, info )
    If ( info < 0) Return

    tst = dds_lapack%state(LAPACK_STATE_Factorized) == MPH_TRUE
    CHCKASSRT( tst, info )
    If ( info < 0) Return

    !-----------------------------------------------------------------
    ! [1] Initialize local variables
    !-----------------------------------------------------------------

    ! Nullify pointers
    Nullify(A_1)
    Nullify(iPIV)

    ! Strategy related
    sym = MAT%sym

    ! Factors related
    n    =  MAT%n
    ldA  =  MAT%ld
    A_1  => MAT%v

    iPIV => dds_lapack%ipiv

    ! Right-hand-side related
    nrhs =  RHS%n
    ldB  =  RHS%ld
    B    => RHS%v
    
    ! Transpose or not ?
    ! 
    ! = 'N':  A * X = B  (No transpose)
    ! = 'T':  A'* X = B  (Transpose)
    ! = 'C':  A'* X = B  (Conjugate transpose = Transpose)

    TRANS = 'N'

    ! Set uplo

    Select Case (sym)
    Case( DM_SYM_IsSPD, DM_SYM_IsSymmetric )

       If ( MAT%stored == DM_STORED_UPPER ) uplo='U'
       If ( MAT%stored == DM_STORED_LOWER ) uplo='L'
       ! convention, see the warning at top of file for details
       If ( MAT%stored == DM_STORED_ALL   ) uplo='U'

    Case Default

       uplo=' ' ! uplo is not needed 

    End Select

    !-----------------------------------------------------------------
    ! [2] Called the right solve routine
    !-----------------------------------------------------------------

    Select Case(sym)
    Case (DM_SYM_IsGeneral)
       
       Call XMPH_ARITHgetrs (& !
            TRANS, N, NRHS,& ! B <-- A^{-1}.B
            A_1, LDA, IPIV,& ! 
            B, LDB, INFO )

    Case (DM_SYM_IsSPD)

       Call XMPH_ARITHpotrs( &
            UPLO, N, NRHS, &
            A_1, LDA,      &
            B, LDB, INFO )

    Case (DM_SYM_IsSymmetric)

       Call XMPH_ARITHsytrs( &
            UPLO, N, NRHS, &
            A_1, LDA, IPIV,&
            B, LDB, INFO )

    Case Default
       info = -1
    End Select

    CHCKASSRT( info >= 0 , info )
    If ( info < 0) Return

  End Subroutine XMPH_LAPACK_solve



End Module XMPH_dds_lapack_mod
