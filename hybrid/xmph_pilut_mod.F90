! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"

#define MPH_DBG_THIS_FILE      0
#define MPH_DBG_FIX_SYMMETRY   1
#define MPH_DBG_IMPORT_MATRIX  0
#define MPH_DBG_IMPORT_PARAM   0
#define MPH_DBG_EXPORT_MATRIX  0
#define MPH_DBG_EXPORT_SMSCHUR 0


!> wrapper to partial ilut, xmph_ilut
Module XMPH_pilut_mod

  ! Modules depencies
  Use Mph_error_mod
  Use XMPH_sparse_matrix_mod

  Implicit None

  !* Derived types *!

  !> structure containing the pilut instance
  Type XMPH_pilut_t
     Sequence

     !> @par Public components
     !------------------------
     !> the input matrix
     Type(XMPH_sparse_matrix_t), Public :: A

     !> the approximation of the schur complement
     Type(XMPH_sparse_matrix_t), Public :: schur
     !> the order of the schur complement
     MPH_INT, Public :: nschur

     !> the limit on LU : number of entries
     MPH_INT, Public :: LUnzRowMax
     !> the limit on the schur complement : number of entries
     MPH_INT, Public :: SnzRowMax

     !> the limit on LU factors : tolerance
     Real(KIND=XMPH_FLOATKIND), Public :: LUtol
     !> the limit on the schur complement : tolerance
     Real(KIND=XMPH_FLOATKIND), Public :: Stol

     !> @par Private components
     !------------------------- 
     !> Internal status
     Integer, Private :: state

  End type XMPH_pilut_t

  !* Defined constants *!

  Integer, Parameter, Private :: HAVE_NOTHING    = 0
  Integer, Parameter, Private :: HAVE_MATRIX     = 1
  Integer, Parameter, Private :: HAVE_SCHURSIZE  = 2
  Integer, Parameter, Private :: HAVE_SCHUR      = 3
  Integer, Parameter, Private :: SCHUR_WAS_GIVEN = 4

  Character(len=MAPHYS_STRL), Private, Parameter :: FLNAME = &
       "XMPH_ARITHmph_pilut_mod.F90"

  !* Defined routines *! 

  Public :: XMPH_pilut_init
  Public :: XMPH_pilut_exit
  Public :: XMPH_pilut_nullify

  Public :: XMPH_pilut_set_matrix
  Public :: XMPH_pilut_set_permutation
  Public :: XMPH_pilut_set_ILULimits
  Public :: XMPH_pilut_set_schurOrder
  Public :: XMPH_pilut_get_schur

  Public :: XMPH_pilut_factorize

  Contains

    ! Constructors and destructors
    !=============================

    !> Init the instance
    !!
    !! @param [inout] this   the pilut instance
    !!
    Subroutine XMPH_pilut_init(this)

      Implicit None
      Type(XMPH_pilut_t), Intent(inout) :: this

      Call XMPH_pilut_nullify(this)

      this%nschur  = 0
      this%LUnzRowMax = 0
      this%SnzRowMax  = 0

      this%LUtol   = 0._XMPH_FLOATKIND
      this%Stol    = 0._XMPH_FLOATKIND

      this%state = HAVE_NOTHING

    End Subroutine XMPH_pilut_init
    !
    !---------------------------------------------------------------------------
    !
    !> Nullify the instance
    !!
    !! nullify every pointers referenced by the instance
    !!
    !! @param [inout] this   the pilut instance
    !!
    Subroutine XMPH_pilut_nullify(this)

      Implicit None
      
      Type(XMPH_pilut_t), Intent(inout) :: this
      Integer :: MPH_INFO_IGNORE

      Call XMPH_sm_nullify(this%A,MPH_INFO_IGNORE)
      Call XMPH_sm_nullify(this%schur,MPH_INFO_IGNORE)

    End Subroutine XMPH_pilut_nullify
    !
    !---------------------------------------------------------------------------
    !
    !> Destroy the instance
    !! 
    !! Destroy the instance by freeing every allocate data,
    !! except the schur complement if it was previous given 
    !! to the user by calling "XMPH_pilut_get_schur)"
    !!
    !! @param [inout] this   the pilut instance
    !!
    Subroutine XMPH_pilut_exit(this)

      Implicit None
      Type(XMPH_pilut_t), Intent(inout) :: this
      Integer :: MPH_INFO_IGNORE

      If (this%state == HAVE_NOTHING) Return

      If (this%state >= HAVE_MATRIX ) &
           Call XMPH_sm_free(this%A, MPH_INFO_IGNORE)

      If (this%state == SCHUR_WAS_GIVEN) Then
         Call XMPH_sm_nullify(this%schur,MPH_INFO_IGNORE)
      Else
         Call XMPH_sm_free(this%schur,MPH_INFO_IGNORE)
      End If

      this%state = HAVE_NOTHING
    End Subroutine XMPH_pilut_exit
    !

    ! Getters and setters
    !====================
    !
    !> Set the matrix of the system
    !! 
    !! Copy a sparse matrix for subsequent factorisation.
    !!
    !! @param [inout] this   the pilut instance
    !! @param [in]    sm     the sparse matrix 
    !! @param [out]   info   the return status
    !!
    Subroutine XMPH_pilut_set_matrix(this,sm,info)

      Implicit None
      Type(XMPH_pilut_t)         , Intent(inout) :: this
      Type(XMPH_sparse_matrix_t) , Intent(in)    :: sm
      Integer                    , Intent(out)   :: info

      Call XMPH_sm_dup(this%A,sm,info)
      CHCKASSRT(info >= 0, info)

    End Subroutine XMPH_pilut_set_matrix
    !
    !---------------------------------------------------------------------------
    !
    !> Set a symmetric permutation on matrix of the system
    !! 
    !!
    !! @param [inout] this   the pilut instance
    !! @param [in]    perm   the permutation
    !! @param [in]    nperm  the permutation size
    !! @param [out]   info   the return status
    !!
    Subroutine XMPH_pilut_set_permutation(this,perm,info)
      Use XMPH_sparse_matrix_mod
      Implicit None
      Type(XMPH_pilut_t)         , Intent(inout) :: this
      MPH_INT, Pointer        , Intent(in)    :: perm(:)
      Integer                    , Intent(out)   :: info

      Call XMPH_sm_convert(this%A,SM_FMT_IJV,info)
      MPH_ONFAILURE_RETURN(info)

      Call XMPH_invperm(this%A%i,perm,this%A%nnz, info)
      MPH_ONFAILURE_RETURN(info)
      
      Call XMPH_invperm(this%A%j,perm,this%A%nnz, info)
      MPH_ONFAILURE_RETURN(info)

      Call XMPH_sm_convert(this%A,SM_FMT_CSR,info)
      MPH_ONFAILURE_RETURN(info)

    End Subroutine XMPH_pilut_set_permutation
    !
    !---------------------------------------------------------------------------
    !
    !> Set the ILU limits 
    !! 
    !! Set the differents limits of the Incomplete factorisation.
    !! Here "S" denotes the schur complement, 
    !! "LU" the L or U part of the factorisation,
    !! "k" the entry index. (since matrices are sparse, in coordinate formats)
    !!
    !! On those objects "O" (=S,L or U), we only store their entries
    !! if they verify 2 assertions "|O(k)| >= tol" and "k <= nzRowMax".
    !!
    !! @param [inout] this     the pilut instance
    !! @param [in]    LUnzRowMax  the maximal number of entries per row in the factors
    !! @param [in]    LUtol    the tolerance in the factors
    !! @param [in]    SnzRowMax   the maximum number of entries per row in the schur
    !! @param [in]    Stol     the tolerance in the schur complement.
    !! @param [out]   info     the return status
    !!
    Subroutine XMPH_pilut_set_ILULimits&
         (this,LUnzRowMax,LUtol,SnzRowMax,Stol,info)

      Implicit None

      !* Arguments *!

      Type(XMPH_pilut_t), Intent(inout)     :: this
      MPH_INT, Intent(in)                :: LUnzRowMax
      Real(KIND=XMPH_FLOATKIND), Intent(in) :: LUtol
      MPH_INT, Intent(in)                :: SnzRowMax
      Real(KIND=XMPH_FLOATKIND), Intent(in) :: Stol
      Integer, Intent(out)                  :: info

      !* Local variables *!

      Integer :: checks(2)
      
      !- End of header ---------------------------------------------------------

      this%LUnzRowMax = LUnzRowMax
      this%SnzRowMax  = SnzRowMax
      this%LUtol   = LUtol
      this%Stol    = Stol

      checks = 0
      CHCKASSRT(LUnzRowMax > 0, checks(1))
      CHCKASSRT(SnzRowMax > 0, checks(2))
      info = Min(checks(1),checks(2))

    End Subroutine XMPH_pilut_set_ILULimits
    !
    !---------------------------------------------------------------------------
    !
    !> Set the Schur complement order
    !! 
    !! Set the schur complement order "schurorder".
    !! Assumes that the schur complement are in the "schurorder" last 
    !! rows/columns of the matrix "this%A"
    !!
    !! @param [inout] this       the pilut instance
    !! @param [in]    schurorder the order of the schur complement
    !! @param [out]   info       the return status
    !!
    Subroutine XMPH_pilut_set_schurorder &
         (this,schurorder,info)
      
      Implicit None
      
      !* Arguments *!

      Type(XMPH_pilut_t) , Intent(inout) :: this
      MPH_INT         , Intent(in   ) :: schurorder
      Integer            , Intent(out  ) :: info
      
      !- End of header ---------------------------------------------------------

      this%nschur = schurorder
      this%state  = HAVE_SCHURSIZE
      info = 0
      CHCKASSRT(schurorder >= 0, info)

    End Subroutine XMPH_pilut_set_schurorder
    !
    !---------------------------------------------------------------------------
    !
    !> Get the Schur complement approximation
    !! 
    !! Get the schur complement matrix.
    !! It must be called after the factorization.
    !!
    !! @param [inout] this       the pilut instance
    !! @param [in]    apschur    the approximation of the schur complement. 
    !!    warning apschur is in csr format.
    !!    and apschur%nnz is not the number of entries
    !!    but the number of allocated elements.
    !! @param [out]   info       the return status
    !!
    Subroutine XMPH_pilut_get_schur &
      (this,apschur,info)
      
      Implicit None

      Type(XMPH_pilut_t)         , Intent(inout) :: this
      Type(XMPH_sparse_matrix_t) , Intent(out)   :: apschur
      Integer                    , Intent(out)   :: info

      info = 0
      CHCKASSRT(this%state >= HAVE_SCHUR, info)
      If (info /= 0) Return

      apschur = this%schur

      this%state = SCHUR_WAS_GIVEN

    End Subroutine XMPH_pilut_get_schur
    !
    ! Methods
    !========
    !> Perform the factorization 
    !! 
    !! Factorize the system and compute the schur complement matrix.
    !! The matrix, the schur order and the pilut parameters 
    !! must be previously set.
    !! Note : Azzam Haidar modified pilut in order to not store the factorization.
    !!
    !! @param [inout] this       the pilut instance
    !! @param [out]   info       the return status
    !!
    Subroutine XMPH_pilut_factorize(this,info)
#if MPH_DBG_THIS_FILE
      Use mph_dbg_mod
#endif
      Use XMPH_sparse_matrix_mod

      Use mph_log_mod
      Implicit None

      !* Arguments *!

      Type(XMPH_pilut_t), Intent(inout) :: this
      Integer, Intent(out) :: info

      !* Local variables *! 

      MPH_INT, Parameter :: UNUSED_INT     = -1
      MPH_INT, Parameter :: UNUSED_PTRSIZE = 1
      Integer  :: MPH_INFO_IGNORED
      Integer  :: itmp

      !* Additional packages *!
      Include 'mpif.h'

      !* Arguments to the routine (by order) *!
      
      MPH_INT :: n
      XMPH_FLOAT, Pointer :: a(:)
      MPH_INT, Pointer :: ja(:)
      MPH_INT, Pointer :: ia(:)

      MPH_INT :: lfil
      Real(KIND=XMPH_FLOATKIND) :: droptol
      MPH_INT :: last

      MPH_INT :: iptr
      MPH_INT :: lalu
      XMPH_FLOAT, Pointer :: alu(:)
      MPH_INT, Pointer :: jlu(:)
      MPH_INT, Pointer :: ju(:)

      XMPH_FLOAT, Pointer :: w(:)
      MPH_INT, Pointer :: jw(:)

      Integer :: ierr

      MPH_INT :: lfilschur
      Real(KIND=XMPH_FLOATKIND) :: dropschur
      Real   (KIND=8) :: timing
      Character(len=MAPHYS_STRL) :: file

      !- End of header ---------------------------------------------------------

      !-------------------------------------------------------------------------
      ! Init 
      !-------------------------------------------------------------------------
#if MPH_DBG_THIS_FILE
      Call MPH_dbg_init
#endif
      Nullify(ia,ja,jlu,ju,jw) ! pointers to integers
      Nullify(a,alu,w)              ! pointers to floats

      !-------------------------------------------------------------------------
      ! Check status
      !-------------------------------------------------------------------------

      ! Job already done
      If ( this%state >= HAVE_SCHUR ) Return

      ! Several fields unset
      info = 0                         
      CHCKASSRT(this%state >= HAVE_SCHURSIZE, info )
      If (info < 0) Goto 9999

      !-------------------------------------------------------------------------
      ! Setup the arguments
      !-------------------------------------------------------------------------

      ! Input matrix

#if MPH_DBG_IMPORT_MATRIX
      Call XMPH_sm_free(this%A,ierr)
      Call MPH_dbg_get_prefix(file)
      file=Trim(file)//"smAILUT.mtx"
      Open(UNIT=10,FILE=file,STATUS = 'OLD',ACTION="read")
      Call XMPH_sm_mmread(this%A,10,ierr)
      Close(10)
      ASSRT(ierr == 0)
      Allocate(this%A%csr(this%A%m+1))
      Call XMPH_sm_ind2ptr(this%A%nnz,this%A%i,this%A%m+1,this%A%csr)
#else

#if MPH_DBG_FIX_SYMMETRY
      ! temporary fix on symmetry
      If (this%A%sym /= SM_SYM_IsGeneral)Then
         Call Mph_log(MSG_WARNING,&
              "Unsymmetrize local matrix.&
              & Reason : PILUT currently only support general matrices")
         Call XMPH_sm_mirror(this%A,info)
         MPH_ONFAILURE_GOTO9999(info)
         this%A%sym = SM_SYM_IsGeneral
         
      End If
#endif
      Call XMPH_sm_convert(this%A, SM_FMT_CSR, info )
      CHCKASSRT(info >= 0, info)
      If (info < 0) Goto 9999
#endif

      n    = this%A%n
      a    => this%A%v
      ja   => this%A%j
      ia   => this%A%csr

#if MPH_DBG_EXPORT_MATRIX
      Call MPH_dbg_init
      Call MPH_dbg_set_file("thisAcsr")
      Call XMPH_sm_mmwrite(this%A,dbg_unit,info) 
#endif

      ! how to drop 

      lfil = this%LUnzRowMax
      droptol = this%LUtol
      last = n - this%nschur

#if MPH_DBG_IMPORT_PARAM
      Call MPH_dbg_get_prefix(file)
      file=Trim(file)//"paramILUT.data"
      Open(UNIT=11,FILE=file,STATUS = 'OLD',ACTION="read")
      Read(11,*) lfil
      Read(11,*) droptol
      Read(11,*) last
      Read(11,*) iptr
      Read(11,*) lalu
      Read(11,*) lfilschur
      Read(11,*) dropschur      
      ! Read(11,*) fullU_entries
      Close(11)
#endif

      ! schur matrix

#if (MPH_DBG_IMPORT_PARAM != 1)
      iptr = 1
      lalu = (2*lfil + 1)*n
#endif

      Call XMPH_sm_create(this%schur, lalu , info )
      CHCKASSRT(info >= 0, info)
      If (info < 0) Goto 9999
      
      this%schur%m   = this%nschur
      this%schur%n   = this%nschur
      this%schur%fmt = SM_FMT_CSR
      Allocate(this%schur%csr(Max(2*n+1,this%schur%n+1)),STAT=info)
      CHCKASSRT(info == 0, info)
      If (info < 0) Goto 9999

      alu => this%schur%v
      jlu => this%schur%j
      ju  => this%schur%csr

      ! workspace
      Allocate( jw(2*n), w(n), STAT=info) 
      CHCKASSRT(info == 0, info)
      If (info < 0) Goto 9999

      ! return status
      ierr = 0
      
      ! new arguments
#if (MPH_DBG_IMPORT_PARAM != 1)
      lfilschur = this%SnzRowMax
      dropschur = this%Stol
#endif

      !-------------------------------------------------------------------------
      ! Call pilut itself
      !-------------------------------------------------------------------------


      Call XMPH_pilut( &
           n,a(1),ja(1),ia(1), lfil,droptol,last, lalu,alu(1),jlu(1),ju(1), &
           w(1),jw(1),ierr, lfilschur,dropschur,timing)
      CHCKASSRT( ierr == 0, info )
      If (info < 0) Call Mph_logWithInfo &
           (MSG_ERROR,ierr,"PILUT, exited with error code:")
      If (info < 0) Goto 9999

      Call XMPH_sm_ptr2ind(this%schur%n+1,ju,this%schur%i)
      this%schur%fmt = SM_FMT_CSR
      this%schur%sym = this%A%sym
      this%schur%nnz = this%schur%csr(this%schur%m+1) - 1

#if (MPH_DBG_EXPORT_SMSCHUR)
      Select Case(this%schur%sym)
      Case(SM_SYM_IsGeneral  ); Call MPH_dbg_set_file("smS_gen.mtx");
      Case(SM_SYM_IsSPD      ); Call MPH_dbg_set_file("smS_spd.mtx");
      Case(SM_SYM_IsSymmetric); Call MPH_dbg_set_file("smS_sym.mtx");
      End Select
      Call XMPH_sm_mmwrite(this%schur, dbg_unit, info)
#endif

#if (MPH_DBG_THIS_FILE)
      Call MPH_dbg_exit
#endif
      !-------------------------------------------------------------------------
      ! Handle errors & free temporary variables 
      !-------------------------------------------------------------------------

9999  Continue

      If ( info == 0 ) this%state = HAVE_SCHUR
      If ( info <  0 ) Call XMPH_sm_free(this%schur, MPH_INFO_IGNORED )
      
      ! Nullify pointers to freed structures & inputs
      Nullify(ia,ja,jlu,ju)
      Nullify(a,alu)
      
      ! Free the rest
      If (Associated(jw )) Deallocate(jw  , STAT=MPH_INFO_IGNORED )
      If (Associated(w  )) Deallocate(w   , STAT=MPH_INFO_IGNORED )

    End Subroutine XMPH_pilut_factorize



End Module XMPH_pilut_mod
