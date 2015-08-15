! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"

#if HAVE_LIBPASTIX
#include "pastix_fortran.h"
#endif

#ifndef MAPHYS_PASTIX_HAS_SETSCHURARRAY
#define MAPHYS_PASTIX_HAS_SETSCHURARRAY 1
#endif

#ifdef MAPHYS_DEBUG
#define MAPHYS_PASTIX_VERBOSE 1
#endif

#ifndef MAPHYS_PASTIX_VERBOSE
#define MAPHYS_PASTIX_VERBOSE 0
#endif

! forgotten index in fortran api
#ifndef IPARM_SCHUR
#define IPARM_SCHUR 28
#endif


! [+] module : sds_pastix_mod --------------------------------------------------
!
!>    MaPHyS' PaSTiX wrapper module
!! @see direct_solver_mod for details
Module XMPH_sds_pastix_mod

#if HAVE_LIBPASTIX
  !* Modules *!
  Use XMPH_sparse_matrix_type

  !* No implicit typing *!
  Implicit None

  !* Type definitions *!

  !> maphys structure for a pastix system
  Type XMPH_pastix_t; Sequence

     !* pastix usual arguments *!
     pastix_data_ptr_t :: data           !< adress of pastix instance
     Integer           :: comm           !< communicator

     Type(XMPH_sparse_matrix_t), pointer :: M !< copy of the input matrix (in CSC)
     pastix_int_t, pointer :: perm(:)    !< permutation tabular
     pastix_int_t, pointer :: invp(:)    !< reverse permutation tabular

     integer               :: nrhs       !< number of right and side (must be 1)
     XMPH_FLOAT, pointer :: rhs (:)    !< right hand side
     XMPH_FLOAT :: dummy_rhs (1)       !< dummy argument for pastix_fortran

     pastix_int_t :: iparm(IPARM_SIZE)   !< Integer parameters
     real(kind=8) :: dparm(DPARM_SIZE)   !< Floating point parameters

     !* thread related *!
     integer , pointer :: bindtab (:)    !< thread binding (size nbthread)

     !* Schur related *!
     Integer                :: nschur    !< size of schur
     Integer                :: ldschur   !< leading dimension
     XMPH_FLOAT, pointer :: schur (:)  !< schur values in dense COLUMN MAJOR format

  End Type XMPH_pastix_t

  !* Private constants *!
  Character(len=MAPHYS_STRL), Parameter, Private :: FLNAME = &
       "XMPH_ARITHmph_sds_pastix_mod.F90"

  !* explicit interfaces *!

  Interface

     !> pastix solver subroutine
     Subroutine XMPH_ARITH_pastix_fortran             &
          (pastix_data,fortran_comm,         &
          n,colptr,row,avals,               &
          perm,invp,b,rhs,iparm,dparm)
       implicit none
       pastix_data_ptr_t:: pastix_data  !# pastix structure's adress
       integer          :: fortran_comm !# fortran communicator

       pastix_int_t     :: n             !# matrix's number of degree of liberty
       pastix_int_t     :: colptr (*)    !# matrix's column in compressed format (n+1)
       pastix_int_t     :: row    (*)    !# matrix's row (colptr(n+1) - 1)
       XMPH_FLOAT     :: avals  (*)    !# matrix's values (colptr(n+1) - 1)

       pastix_int_t     :: perm   (*)    !# permutations 
       pastix_int_t     :: invp   (*)    !# inverse permutations

       XMPH_FLOAT     :: b      (*)    !# right-hand sides values (rhs * n) 
       pastix_int_t     :: rhs           !# number of right hand side            
       integer          :: iparm  (*)    !# pastix's integer parameters
       real(kind=8)     :: dparm  (*)    !# pastix's float   parameters
     End Subroutine XMPH_ARITH_pastix_fortran

     !> check the csc matrix and eventually ask to reallocate the matrix (symetrise it for scotch)
     Subroutine XMPH_ARITH_pastix_fortran_checkmatrix &
          (data_check, fortran_comm, verb, flagsym, flagcor, &
          n, colptr, row, avals, loc2glob, dof)
       implicit none
       pastix_data_ptr_t :: data_check    !# pastix check data structure adress   
       integer           :: fortran_comm  !# fortran communicator 
       integer           :: verb          !# pastix verbosity (6: maximum) 

       integer           :: flagsym       !# matrix' symmetry  
       integer           :: flagcor       !# indicate if we permit the function 
       !                                  !# to reallocate the matrix. 

       pastix_int_t      :: n             !# size of the matrix
       pastix_int_t      :: colptr (*)    !# matrix's column in 
       !# compressed format (n+1)
       pastix_int_t      :: row    (*)    !# matrix's row (colptr(n+1) - 1)
       XMPH_FLOAT      :: avals  (*)    !# matrix's values (colptr(n+1) - 1)

       pastix_int_t      :: loc2glob      !# global column number of 
       !                                  !# the local columns,
       !                                  !#  -1 if not distributed.
       integer           :: dof           !# number of degrees of freedom.

     End Subroutine XMPH_ARITH_pastix_fortran_checkmatrix

     !> symetrise the matrix (for scotch)
     subroutine XMPH_ARITH_pastix_fortran_checkmatrix_end &
          (data_check, verb, row, avals, dof)
       implicit none
       pastix_data_ptr_t :: data_check    !# pastix check data structure adress  
       integer           :: verb          !# pastix verbosity (6: maximum) 
       pastix_int_t      :: row    (*)    !# matrix's row     (size colptr(n+1) - 1)
       XMPH_FLOAT      :: avals  (*)    !# matrix's values  (size colptr(n+1) - 1)
       integer           :: dof           !# number of degrees of freedom.

     End Subroutine XMPH_ARITH_pastix_fortran_checkmatrix_end

  End Interface

#else
  ! Module to log error messages
  Use mph_log_mod
  Implicit None

  ! Dummy type
  Type XMPH_pastix_t; Sequence
     Integer :: dummy
  End type XMPH_pastix_t

  ! error constants
  Integer, Parameter, Private :: ERR_LIBPASTIX_UNAVAILABLE = -9999
  Character(len=MAPHYS_STRL), Parameter, Private :: &
       LIBPASTIX_UNAVAILABLE = "library libpastix is unavailable"
#endif

  !* Access specifiers *!
  Public :: XMPH_pastix_available
  Public ::  XMPH_pastix_set_MPIcommunicator 
  Public ::  XMPH_pastix_set_matrix 
  Public ::  XMPH_pastix_analyze 
  Public ::  XMPH_pastix_factorize 
  Public ::  XMPH_pastix_solve_RHS  
  Public ::  XMPH_pastix_finalize 
  Public ::  XMPH_pastix_set_denseschur  
  Public ::  XMPH_pastix_solve_sparseRHS 
  Public ::  XMPH_pastix_set_distributed_matrix 
  Public ::  XMPH_pastix_set_multithread 
  Public ::  XMPH_pastix_get_permutation
  Public ::  XMPH_pastix_get_numstats
  Public ::  XMPH_pastix_get_memstats
  Public ::  XMPH_pastix_get_perfstats
  Public ::  XMPH_pastix_set_schurlist
  Public ::  XMPH_pastix_get_schur
  Public ::  XMPH_pastix_free_schur

  Private ::  XMPH_pastix_exec
  Private ::  XMPH_pastix_check_status
  Private ::  XMPH_pastix_sparse_matrix_check
#if HAVE_LIBPASTIX
  Private ::  XMPH_ARITH_pastix_fortran       
  Private ::  XMPH_ARITH_pastix_fortran_checkmatrix 
  Private ::  XMPH_ARITH_pastix_fortran_checkmatrix_end
#endif



  !* routines *!

Contains

  ! [+] function : XMPH_pastix_available ---------------------------------------
  ! 
  !> Return 0 if pastix library is available or -1 if not.
  Function XMPH_pastix_available()
    Integer :: XMPH_pastix_available
#if HAVE_LIBPASTIX    
    XMPH_pastix_available = 0
#else
    XMPH_pastix_available = -1
#endif
  End Function XMPH_pastix_available



  ! [+] routine : XMPH_pastix_set_MPIcommunicator ------------------------------
  !
  !> Set the MPI communicator to be used
  !! @see XMPH_sds_set_MPIcommunicator()
  Subroutine XMPH_pastix_set_MPIcommunicator (id_pastix, comm, info  )

    !* Modules & co *!
    Use mph_error_mod
    Use mph_log_mod
    Implicit None
    Include 'mpif.h'

    !* Arguments *!
    Type(XMPH_pastix_t), intent(inout) :: id_pastix
    Integer              , intent(in   ) :: comm
    Integer              , intent(  out) :: info 

#if HAVE_LIBPASTIX    

    !* Local variables *!
    ! Constants
    Integer, Parameter :: INI_VAL=0
    Integer, Parameter :: INT_EQUIVALENT_TO_NULL_PTR = 0

    ! Scalars
    Integer               :: iinfo

    !- End of header -----------------------------------------------------------

    info =0

    !---------------------------------------------------------------------------
    ! [1] Init id_pastix structure
    !---------------------------------------------------------------------------

    !> @warning :
    !! id_pastix%data must contain the adress to the NULL pointer.
    !! Indeed Pastix schematically do : 
    !! if ( id_pastix%data ==  NULL) Initialize_pastix(id_pastix%data)
    !!
    id_pastix%data = INT_EQUIVALENT_TO_NULL_PTR ! set to NULL
    id_pastix%comm = INI_VAL
    id_pastix%nrhs  = INI_VAL
    id_pastix%nschur  = INI_VAL
    id_pastix%ldschur = INI_VAL

    id_pastix%iparm = INI_VAL
    id_pastix%dparm = INI_VAL
    id_pastix%dummy_rhs = INI_VAL
    !
    Nullify(id_pastix%M)
    Nullify(id_pastix%perm)
    Nullify(id_pastix%invp)
    Nullify(id_pastix%rhs)
    Nullify(id_pastix%bindtab)
    Nullify(id_pastix%schur)

    !---------------------------------------------------------------------------
    ! [2] set up the MPI communicator
    !---------------------------------------------------------------------------

    Call MPI_Comm_dup(comm, id_pastix%comm , iinfo )
    ASSRT(iinfo == MPI_SUCCESS)

#else
    info = ERR_LIBPASTIX_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBPASTIX_UNAVAILABLE )
    Return
#endif 


  End Subroutine XMPH_pastix_set_MPIcommunicator

  ! [+] routine : XMPH_pastix_set_matrix ---------------------------------------
  !
  !> Set the matrix centralized on process 0.
  !! @see XMPH_sds_set_matrix()
  Subroutine XMPH_pastix_set_matrix (id_pastix, M, info )

    !* Modules *!
    Use XMPH_sparse_matrix_mod
    Use mph_error_mod
    Use mph_log_mod
    Implicit None
    Include 'mpif.h'

    !* Arguments *!

    Type(XMPH_pastix_t), Intent(inout) :: id_pastix
    Type(XMPH_sparse_matrix_t), Intent(in   ) :: M
    Integer              , Intent(  out) :: info

#if HAVE_LIBPASTIX
    !* Local variables *!

    ! Scalars

    Integer :: rank    
    Integer :: COMM_SELF ! copy of MPI_COMM_SELF (necessary)
    MPH_INT :: i              
    Integer :: verbose
    Logical :: maySegfault

    ! dummy
    Type(XMPH_sparse_matrix_t), pointer:: A           !< pointer to id_pastix%M

    !-End of header ------------------------------------------------------------

    !> @note COMM_SELF. Why do you need it ?
    !! In XMPH_pastix_sparse_matrix_check()
    !! the communicator is an "intent(inout)"
    !! .

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------

    info=0

    Call MPI_Comm_rank(id_pastix%comm, rank, info) 
    ASSRT(info == MPI_SUCCESS)

    Call mph_log_GetVerbosity(verbose)
    If (verbose >= MSG_DEBUG ) Then
       verbose = 6
    Else
       verbose = API_VERBOSE_NOT
    End If

    ! perform some checking which may lead to segfault.
    ! If pastix runs on near-dense matrix and the
    ! number of elements is huge, segfault may occurs.
    maySegfault =( M%m*M%n /= M%nnz ).And.(M%nnz >= Huge(M%nnz)/(MPH_INTKIND/2))
    If( maySegfault )&
         Call MPH_LogWithInfo(MSG_WARNING,M%nnz,&
         "(pastix wrapper) matrix is not dense,&
         & and its nnz is near the limit of the integer.&
         & Segmentation fault may occur. nnz =")

    !---------------------------------------------------------------------------
    ! [2] Create the Matrix for PaSTiX
    !---------------------------------------------------------------------------

    ! Pastix Matrix, must be the same on all MPI process,
    ! In CSC format.

    ! allocate the structure to hold the PaSTiX matrix 
    If ( Associated(id_pastix%M) ) Deallocate( id_pastix%M )
    Allocate(id_pastix%M,stat=info);
    CHCKASSRT(info == 0, info)
    If ( info < 0 ) Return

    ! master convert M into pastix format
    master_create_CSC_matrix:If ( rank == 0 ) then

       ! duplicate
       Call XMPH_sm_dup(id_pastix%M,M, info);
       CHCKASSRT(info == 0, info)
       If ( info < 0 ) Return

       ! check the matrix
       COMM_SELF = MPI_COMM_SELF
       call XMPH_pastix_sparse_matrix_check&
            (id_pastix%M,COMM_SELF,info)
       CHCKASSRT(info == 0, info)
       If ( info < 0 ) Return

    End If master_create_CSC_matrix

    !  broadcast the matrix
    Call XMPH_sm_bcast(id_pastix%M,0,id_pastix%comm,info)
    ASSRT(info == MPI_SUCCESS)
    A => id_pastix%M

    !---------------------------------------------------------------------------
    ! [3] initiate Pastix parameters 
    ! [ - ] (only touch id_pastix%iparm and id_pastix%dparm) 
    !---------------------------------------------------------------------------

    !
    Allocate(id_pastix%perm(A%n),stat=info)
    CHCKASSRT(info == 0, info)
    If ( info < 0 ) Return

    Allocate(id_pastix%invp(A%n),stat=info) 
    CHCKASSRT(info == 0, info)
    If ( info < 0 ) Return

    ! Give initial values
    Do i=1, A%n
       id_pastix%perm(i) = 0
       id_pastix%invp(i) = 0
    End Do

    ! Call Pastix to fill iparm/dparm with their default values
    id_pastix%iparm(IPARM_MODIFY_PARAMETER) = API_NO
    Call XMPH_pastix_exec(id_pastix,useRHS=.False.)
    Call XMPH_pastix_check_status(id_pastix, __LINE__, info)
    If ( info < 0 ) Return

    !---------------------------------------------------------------------------
    ! [4] Customize parameters iparm/dparm
    !---------------------------------------------------------------------------

    id_pastix%iparm(IPARM_VERBOSE)          = verbose
    id_pastix%iparm(IPARM_START_TASK)       = API_TASK_INIT
    id_pastix%iparm(IPARM_END_TASK)         = API_TASK_INIT

    id_pastix%dparm(DPARM_EPSILON_MAGN_CTRL)= 1.d-50 ! lty testing pivoting

    Select Case(A%sym)
    Case (SM_SYM_IsSPD)
       id_pastix%iparm(IPARM_SYM)           = API_SYM_YES
       id_pastix%iparm(IPARM_FACTORIZATION) = API_FACT_LLT

    Case (SM_SYM_IsSymmetric)
       id_pastix%iparm(IPARM_SYM)           = API_SYM_YES
       id_pastix%iparm(IPARM_FACTORIZATION) = API_FACT_LDLT

    Case (SM_SYM_IsGeneral )
       id_pastix%iparm(IPARM_SYM)           = API_SYM_NO
       id_pastix%iparm(IPARM_FACTORIZATION) = API_FACT_LU

    End Select

    id_pastix%iparm(IPARM_MATRIX_VERIFICATION) = API_NO
    id_pastix%iparm(IPARM_RHS_MAKING)          = API_RHS_B

    !---------------------------------------------------------------------------
    ! [4] initialise XMPH_pastix_data structure. (id_pastix%data)
    !---------------------------------------------------------------------------

    Call XMPH_pastix_exec(id_pastix,useRHS=.False.)
    Call XMPH_pastix_check_status(id_pastix, __LINE__ , info )
    If (info < 0) Return

    !---------------------------------------------------------------------------
    ! [5] Finish
    !---------------------------------------------------------------------------

    ! deactivate multi threading by default
    id_pastix%iparm(IPARM_THREAD_NBR) = 1
    id_pastix%iparm(IPARM_BINDTHRD)   = API_BIND_NO

    If (maySegfault) info = 1

#else
    info = ERR_LIBPASTIX_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBPASTIX_UNAVAILABLE )
    Return
#endif 


  End Subroutine XMPH_pastix_set_matrix

  ! [+] routine : XMPH_pastix_analyze ------------------------------------------
  ! 
  !> Analyze the linear system.
  !! @see XMPH_sds_analyze()
  subroutine XMPH_pastix_analyze (id_pastix, info  )

    !* Module(s) *!
    Use mph_error_mod
    Implicit None

    !* Argument(s) *!

    Type(XMPH_pastix_t), intent(inout) :: id_pastix
    Integer              , intent(  out) :: info 

#if HAVE_LIBPASTIX
    !* Local variable(s) *!
    
    ! scalar
    MPH_INT :: k

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Call the analysis step
    !---------------------------------------------------------------------------

    id_pastix%iparm(IPARM_ORDERING)   = API_ORDER_SCOTCH
    id_pastix%iparm(IPARM_START_TASK) = API_TASK_ORDERING
    id_pastix%iparm(IPARM_END_TASK)   = API_TASK_ANALYSE   

    Call XMPH_pastix_exec(id_pastix,useRHS=.False.)
    Call XMPH_pastix_check_status(id_pastix, __LINE__, info )
    If (info < 0) Return

    !---------------------------------------------------------------------------
    ! [2] Modify permutation vectors 
    !---------------------------------------------------------------------------
    ! as pastix begin indexes at 0 (C)
    ! and we    begin indexes at 1 (Fortran)

    Do k=1,id_pastix%M%n
       id_pastix%perm(k) = id_pastix%perm(k) + 1
       id_pastix%invp(k) = id_pastix%invp(k) + 1
    End Do

#else
    info = ERR_LIBPASTIX_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBPASTIX_UNAVAILABLE )
    Return
#endif 
    
  End subroutine XMPH_pastix_analyze

  ! [+] routine : XMPH_pastix_factorize ----------------------------------------
  ! 
  !> Factorize the matrix (previously set)
  !! @see XMPH_sds_factorize()
  Subroutine XMPH_pastix_factorize (id_pastix, info  )

    !* Module(s) *!

    Use XMPH_dense_matrix_mod
    Use mph_error_mod
    Use mph_log_mod
    Implicit None

    !* Argument(s) *!
    
    Type(XMPH_pastix_t), intent(inout) :: id_pastix
    Integer              , intent(  out) :: info 
    
#if HAVE_LIBPASTIX
    !* Local variable(s) *!

    ! scalars

    MPH_INT :: nschur, mschur

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [0] Allocate & set local schur complement if needed
    !---------------------------------------------------------------------------


    If ( (id_pastix%iparm(IPARM_SCHUR) == API_YES).And.&
         (.Not. Associated(id_pastix%schur)      ))Then
       ! get the number of columns of the Schur C. on this MPI process
       Call XMPH_ARITH_pastix_fortran_getschurlocalnodenbr( &
            id_pastix%data, nschur, info )
       CHCKASSRT( info == NO_ERR, info )
       If (info < 0) Return
       id_pastix%nschur = nschur

       ! get the number of lines of the Schur C. on this MPI process
       mschur = id_pastix%ldschur

       ! allocate the Schur C.
       Allocate (id_pastix%schur(mschur*nschur),STAT=info)
       CHCKASSRT( info == NO_ERR, info )
       If (info < 0) Return

       ! set it
       Call XMPH_ARITH_pastix_fortran_setschurarray &
            (id_pastix%data, id_pastix%schur(1), info)
       CHCKASSRT( info >= 0, info)
       If ( info < 0) Return

    End If

    !---------------------------------------------------------------------------
    ! [1] Perform factorization 
    !---------------------------------------------------------------------------
    
    id_pastix%iparm(IPARM_START_TASK) = API_TASK_NUMFACT
    id_pastix%iparm(IPARM_END_TASK)   = API_TASK_NUMFACT

    MKL_SET_NUM_THREADS(1) 
    Call XMPH_pastix_exec(id_pastix,useRHS=.False.)
    Call XMPH_pastix_check_status(id_pastix, __LINE__ , info)
    MKL_SET_NUM_THREADS(id_pastix%iparm(IPARM_THREAD_NBR))

    If (info < 0) Return

#else
    info = ERR_LIBPASTIX_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBPASTIX_UNAVAILABLE )
    Return
#endif 

  End Subroutine XMPH_pastix_factorize

  ! [+] routine : XMPH_pastix_solve_RHS ----------------------------------------
  ! 
  !> Solve the linear system with a dense, centralized right-hand-side
  !! @see XMPH_sds_solve_RHS()
  subroutine XMPH_pastix_solve_RHS  &
       (id_pastix, rhs, nrhs, ldrhs, info  )

    !* Module(s) *!
    Use mph_error_mod
    Use mph_log_mod
    Implicit None
    Include 'mpif.h'

    !* Argument(s) *!
    Type(XMPH_pastix_t), intent(inout) :: id_pastix
    Integer              , intent(in   ) :: nrhs,ldrhs
    XMPH_FLOAT     , target   , intent(inout) :: rhs(ldrhs*nrhs)
    Integer              , intent(  out) :: info 
    
#if HAVE_LIBPASTIX
    !* Local variable(s) *!
    
    ! scalars
    Integer :: me, np

    ! type
    Type(XMPH_sparse_matrix_t), pointer:: A !< pointer to id_pastix%M
    
    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Broadcast the RHS
    !---------------------------------------------------------------------------

    A => id_pastix%M

    !---------------------------------------------------------------------------
    ! [1.1] Get MPI data
    !---------------------------------------------------------------------------

    Call MPI_Comm_size(id_pastix%comm, np, info ) 
    ASSRT(info  == MPI_SUCCESS)

    Call MPI_Comm_rank(id_pastix%comm, me, info ) 
    ASSRT(info  == MPI_SUCCESS)

    !---------------------------------------------------------------------------
    ! [1.2] Link/allocate rhs
    !---------------------------------------------------------------------------

    Nullify(id_pastix%rhs)
    id_pastix%nrhs = nrhs
    if (me == 0) id_pastix%rhs => rhs(1:nrhs*A%n)
    if (me /= 0) then

       Allocate(id_pastix%rhs(nrhs*A%n),STAT=info )
       CHCKASSRT( info == 0, info )
       If (info < 0) Goto 9999

    endif

    !---------------------------------------------------------------------------
    ! [1.3] Broadcast RHS
    !---------------------------------------------------------------------------

    Call MPI_Bcast &
         ( id_pastix%rhs(1), nrhs*A%n,&
         & XMPH_FLOATMPI, 0, id_pastix%comm, info )
    ASSRT(info  == MPI_SUCCESS)

    !---------------------------------------------------------------------------
    ! [2] Launch the solve
    !---------------------------------------------------------------------------
        
    ! set the parameters
    id_pastix%iparm(IPARM_START_TASK)       = API_TASK_SOLVE
    id_pastix%iparm(IPARM_END_TASK)         = API_TASK_REFINE  

    ! disallow refinement when computing the schur
    If ( id_pastix%iparm(IPARM_SCHUR) == API_YES ) &
         id_pastix%iparm(IPARM_END_TASK) = API_TASK_SOLVE

    ! call pastix
    MKL_SET_NUM_THREADS(1)
    Call XMPH_pastix_exec(id_pastix,useRHS=.True.)
    Call XMPH_pastix_check_status(id_pastix, __LINE__ , info )
    MKL_SET_NUM_THREADS(id_pastix%iparm(IPARM_THREAD_NBR))

    If ( info < 0) Goto 9999

    !---------------------------------------------------------------------------
    ! [3] Finish
    !---------------------------------------------------------------------------

    info = 0

9999 Continue

    ! Unlink data / Free memory
    id_pastix%nrhs = 0
    If (Associated(id_pastix%rhs))Then
       if (me /= 0)  Deallocate(id_pastix%rhs)
       if (me == 0)  Nullify   (id_pastix%rhs)
    End If
    Nullify (A)

#else
    !---------------------------------------------------------------------------
    ! [*] Libpastix unavailable
    !---------------------------------------------------------------------------
    info = ERR_LIBPASTIX_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBPASTIX_UNAVAILABLE )
    Return
#endif 

  End subroutine XMPH_pastix_solve_RHS

  ! [+] routine : XMPH_pastix_finalize -----------------------------------------
  ! 
  !> Free the pastix instance
  !! @see XMPH_sds_finalize()
  Subroutine XMPH_pastix_finalize (id_pastix, info  )

    !* Modules *!

    Use XMPH_sparse_matrix_mod
    Use mph_error_mod
    Implicit None
    include 'mpif.h'

    !* Arguments *!

    type(XMPH_pastix_t), intent(inout) :: id_pastix
    integer              , intent(  out) :: info 

#if HAVE_LIBPASTIX
    !* Local variables *!
    Type(XMPH_sparse_matrix_t), pointer:: A !< pointer to id_pastix%M

    !- End of header -----------------------------------------------------------

    !--------------------------------------------------------------------------- 
    ! [1] Clean id_pastix%data
    !--------------------------------------------------------------------------- 

    A => id_pastix%M
    id_pastix%iparm(IPARM_START_TASK)       = API_TASK_CLEAN
    id_pastix%iparm(IPARM_END_TASK)         = API_TASK_CLEAN

    !yohan>
    ! id_pastix%rhs must not be null in XMPH_ARITH_pastix_fortran call.
    ! It generates an error with intel compilers 
    !yohan<
    If (id_pastix%nrhs <= 0) Then
       Allocate(id_pastix%rhs(1))
    End If

    Call XMPH_pastix_exec(id_pastix,useRHS=.True.)
    Call XMPH_pastix_check_status(id_pastix, __LINE__, info  )
    If (info  < 0) Return

    !--------------------------------------------------------------------------- 
    ! [2] free memory inside the structure
    !--------------------------------------------------------------------------- 

    Nullify(A)
    Call MPI_Comm_free(id_pastix%comm,info)
    ASSRT(info == MPI_SUCCESS)

    Call XMPH_sm_free(id_pastix%M,info) 
    CHCKASSRT(info >= 0, info)
    If (info < 0) Return

    If(Associated(id_pastix%M   )) Deallocate(id_pastix%M)
    If(Associated(id_pastix%perm)) Deallocate(id_pastix%perm)
    If(Associated(id_pastix%invp)) Deallocate(id_pastix%invp)

    If(Associated(id_pastix%schur)) Deallocate(id_pastix%schur)


    If (id_pastix%nrhs <= 0) Then
       If(Associated(id_pastix%rhs)) Deallocate(id_pastix%rhs)
    Else
       If(Associated(id_pastix%rhs)) Nullify(id_pastix%rhs)
    End If

#else
    !---------------------------------------------------------------------------
    ! [*] Libpastix unavailable
    !---------------------------------------------------------------------------
    info = ERR_LIBPASTIX_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBPASTIX_UNAVAILABLE )
    Return
#endif 

  End Subroutine XMPH_pastix_finalize

  ! [+] routine : XMPH_pastix_set_denseschur -----------------------------------
  ! 
  !> Ask to compute the dense schur during the factorization
  !! @see XMPH_sds_set_denseshur()
  !! 
  !! Warning : see README about specific calls 
  subroutine XMPH_pastix_set_denseschur  &
       (id_pastix, schurlist, nschurlist, schur, nschur, ldschur, info )

    !* Modules *!

    Use mph_error_mod
#if MAPHYS_DEBUG
    Use mph_dbg_mod
#endif
    implicit none

    !* Arguments *!

    Type(XMPH_pastix_t), intent(inout) :: id_pastix
    MPH_INT           , intent(in   ) :: nschurlist
    MPH_INT           , intent(in   ) :: schurlist(nschurlist)
    MPH_INT           , intent(in   ) :: nschur,ldschur
    XMPH_FLOAT, pointer, intent(inout) :: schur(:)
    integer              , intent(  out) :: info 

#if HAVE_LIBPASTIX
    !* Local Variables *!

#if MAPHYS_DEBUG
    MPH_INT :: i
#endif

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] check data
    !---------------------------------------------------------------------------

    CHCKASSRT( nschurlist == nschur, info)
    If ( info < 0) Return
    CHCKASSRT( nschur    == ldschur, info)
    If ( info < 0) Return

    !---------------------------------------------------------------------------
    ! [2] set the schur list 
    !---------------------------------------------------------------------------
#if MAPHYS_DEBUG
    Call MPH_dbg_init

    Call MPH_dbg_set_file("schurlist.mtx")
    Write(dbg_unit,*) nschurlist
    Do i = 1, nschurlist
       Write(dbg_unit,*) schurlist(i)
    End Do

    Call MPH_dbg_set_file("end.txt")
    Call MPH_dbg_barrier
#endif

    Call XMPH_ARITH_pastix_fortran_setschurunknownlist &
         (id_pastix%data, nschurlist, schurlist(1), info )
    CHCKASSRT( info >= 0, info)
    If ( info < 0) Return

    !---------------------------------------------------------------------------
    ! [3] activate the schur option 
    !---------------------------------------------------------------------------

    id_pastix%iparm(IPARM_SCHUR)               = API_YES
    id_pastix%iparm(IPARM_MATRIX_VERIFICATION) = API_NO

    !---------------------------------------------------------------------------
    ! [4] user defined rhs
    !---------------------------------------------------------------------------

    id_pastix%iparm(IPARM_RHS_MAKING)          = API_RHS_B

    !---------------------------------------------------------------------------
    ! [5] set links to the output schur (<XMPH_pastix_factorize>)
    !---------------------------------------------------------------------------

    id_pastix%nschur  =  nschur
    id_pastix%ldschur = ldschur
    id_pastix%schur  =>  schur

#if MAPHYS_PASTIX_HAS_SETSCHURARRAY
    Call XMPH_ARITH_pastix_fortran_setschurarray(id_pastix%data, schur(1), info)
    CHCKASSRT( info >= 0, info)
    If ( info < 0) Return
#endif

#else
    !---------------------------------------------------------------------------
    ! [*] Libpastix unavailable
    !---------------------------------------------------------------------------
    info = ERR_LIBPASTIX_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBPASTIX_UNAVAILABLE )
    Return
#endif 

  End Subroutine XMPH_pastix_set_denseschur


  ! [+] routine : XMPH_pastix_solve_sparseRHS ----------------------------------
  ! 
  !> Solve the linear system with a sparse right-hand-side in coordinate format
  !! @see XMPH_sds_solve_sparseRHS()
  !!
  !! @warning here I optimize the communications and not the number of operations.
  !!
  !! @see direct_solver_solve_sparseRHS
  !!
  !! @author Yohan Lee-tin-yien
  Subroutine XMPH_pastix_solve_sparseRHS &
       (id_pastix, nx, ldx, x, sp_rhs, info )

    !* Module(s) *!

    Use XMPH_sparse_matrix_mod
    Use mph_error_mod
    Implicit None
    Include 'mpif.h'
    
    !* Argument(s) *!

    Type(XMPH_pastix_t), intent(inout) :: id_pastix
    Integer              , intent(in   ) :: nx,ldx
    XMPH_FLOAT , target, intent(inout) :: x(nx*ldx)
    Type(XMPH_sparse_matrix_t), intent(in   ) :: sp_rhs ! in CSC format
    Integer              , intent(  out) :: info 

#if HAVE_LIBPASTIX
    !* Local variable(s) *!

    Integer                       :: me
    Type(XMPH_sparse_matrix_t), pointer:: A !< pointer to id_pastix%M
    Type(XMPH_sparse_matrix_t)         :: sp_rhs_cpy !< copy of sp_rhs
    XMPH_FLOAT, pointer              :: rhs1 (:) ! temporary rhs
    integer                       :: i, j, k

    !- End of header -----------------------------------------------------------
    
    !---------------------------------------------------------------------------
    ! [1] Check arguments
    !---------------------------------------------------------------------------
    
    CHCKASSRT(ldx == sp_rhs%m, info)
    If ( info < 0 ) Return

    CHCKASSRT( nx == sp_rhs%n, info)
    If ( info < 0 ) Return

    !---------------------------------------------------------------------------
    ! [2] Broadcast the sparse RHS
    !---------------------------------------------------------------------------

    Call MPI_Comm_rank( id_pastix%comm, me, info ) 
    ASSRT( info == MPI_SUCCESS)

    If (me == 0)Then
       Call XMPH_sm_dup(sp_rhs_cpy,sp_rhs,info )
       CHCKASSRT( info >= 0, info)
       If ( info < 0 ) Return
    End If

    Call XMPH_sm_bcast(sp_rhs_cpy, 0, id_pastix%comm, info )
    ASSRT( info == MPI_SUCCESS)

    !---------------------------------------------------------------------------
    ! [3] Convert sparse RHS into a dense RHS
    !---------------------------------------------------------------------------

    id_pastix%nrhs = nx
    Nullify(id_pastix%rhs)
    Allocate(id_pastix%rhs(sp_rhs_cpy%m * sp_rhs_cpy%n), STAT=info )
    CHCKASSRT( info == 0, info)
    If ( info < 0 ) Return
    
    Do k= 1, sp_rhs_cpy%m * sp_rhs_cpy%n
       id_pastix%rhs(k) = XMPH_FLOATZERO
    End Do

    do k=1,sp_rhs_cpy%nnz
       i = sp_rhs_cpy%i(k)
       j = sp_rhs_cpy%j(k)
       id_pastix%rhs( i + (j-1)*sp_rhs_cpy%m ) = sp_rhs_cpy%v(k) 
    enddo
    Call XMPH_sm_free(sp_rhs_cpy,info )
    CHCKASSRT( info >= 0, info)
    If ( info < 0 ) Return

    !---------------------------------------------------------------------------
    ! [4] Launch the Solve
    !---------------------------------------------------------------------------

    Do k= 1, nx*ldx
       x(k) = XMPH_FLOATZERO
    End Do

    A => id_pastix%M
    Allocate(rhs1(ldx), STAT= info )
    CHCKASSRT( info == 0, info)
    If ( info < 0 ) Return

    do j=1, nx ! for each right hand side, perform a solve+refine

       Do k=1,ldx
          rhs1(k) = id_pastix%rhs( k + (j-1)*ldx)
       End Do

       id_pastix%iparm(IPARM_START_TASK)       = API_TASK_SOLVE
       id_pastix%iparm(IPARM_END_TASK)         = API_TASK_REFINE  

       MKL_SET_NUM_THREADS(1)
       Call XMPH_ARITH_pastix_fortran              &
            (id_pastix%data, id_pastix%comm, &
            A%n,A%csc,A%i,A%v,              &
            id_pastix%perm,id_pastix%invp,   &
            rhs1        ,    1        ,    &
            id_pastix%iparm,id_pastix%dparm )
       Call XMPH_pastix_check_status(id_pastix, __LINE__, info )
       MKL_SET_NUM_THREADS(id_pastix%iparm(IPARM_THREAD_NBR))

       If ( info < 0 ) Return

       ! save the solution
       If (me == 0 )Then
          Do k=1,ldx
             x( k + (j-1)*ldx ) = rhs1(k)
          End Do
       End If
    End do

    !---------------------------------------------------------------------------
    ! [5] Free memory
    !---------------------------------------------------------------------------

    ! free memory   
    nullify(A)
    deallocate(rhs1)
    deallocate(id_pastix%rhs)

#else
    !---------------------------------------------------------------------------
    ! [*] Libpastix unavailable
    !---------------------------------------------------------------------------
    info = ERR_LIBPASTIX_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBPASTIX_UNAVAILABLE )
    Return
#endif 

  End subroutine XMPH_pastix_solve_sparseRHS

  ! [+] routine : XMPH_pastix_set_distributed_matrix ---------------------------
  !
  !> set pastix instance with a distributed matrix in coordinate format.
  !!
  !! @param id_pastix pastix structure
  !! @param M         the distributed matrix 
  !! @param global_n  the order of the distributed matrix
  !! @param info      the status of the subroutine
  !!
  !! @note feature do not exist in pastix.
  !! pastix does not support distributed matrix in ijv format
  !! so matrix M is "Allgathered" and assembled into matrix Mnew.
  !! Finally, XMPH_pastix_set_matrix is called with Mnew .
  !!
  subroutine XMPH_pastix_set_distributed_matrix &
       (id_pastix, M, global_n, info  )

    !* Module(s) *!

    Use XMPH_sparse_matrix_mod
    Use mph_error_mod
#if MAPHYS_DEBUG
    Use mph_dbg_mod
#endif
    Implicit none
    Include 'mpif.h'

    !* Argument(s) *!

    type(XMPH_pastix_t), intent(inout) ::  id_pastix
    type(XMPH_sparse_matrix_t), intent(in   ) ::  M       
    integer              , intent(in   ) :: global_n 
    integer              , intent(  out) :: info     
    
#if HAVE_LIBPASTIX
    !* Local variable(s) *!
    
    ! constants
    Integer, Parameter :: root = 0
#if MAPHYS_XMPH_PASTIX_VERBOSE
    Integer, Parameter :: verbose = API_VERBOSE_YES
#else
    Integer, Parameter :: verbose = API_VERBOSE_NO
#endif

    ! types
    Type(XMPH_sparse_matrix_t), Pointer:: A           !< pointer to id_pastix%M

    !- End of header -----------------------------------------------------------

    !--------------------------------------------------------------------------- 
    ! [1] All gather the matrix
    !--------------------------------------------------------------------------- 

    Allocate(id_pastix%M, STAT = info )
    CHCKASSRT(info == 0, info)
    If (info<0) Return

    Call XMPH_sm_Allgather &
         (M, id_pastix%M, global_n, id_pastix%comm, info )
    CHCKASSRT(info == 0, info)
    If (info<0) Return


#ifdef MAPHYS_DEBUG
    Call MPH_dbg_init
    Call MPH_dbg_set_file("sm_new.mtx")
    Call XMPH_sm_mmwrite(id_pastix%M, dbg_unit, info )
#endif

    !--------------------------------------------------------------------------- 
    ! [2] Convert to CSC
    !--------------------------------------------------------------------------- 

    Call XMPH_sm_assemble &
         (id_pastix%M, global_n, info) 
    CHCKASSRT(info == 0, info)
    If (info<0) Return

    Call XMPH_sm_convert(id_pastix%M, SM_FMT_CSC,info)
    CHCKASSRT(info >= 0, info)
    If (info<0) Return

    Call XMPH_pastix_sparse_matrix_check &
         (id_pastix%M,id_pastix%comm, info)
    CHCKASSRT(info >= 0, info)
    If (info<0) Return

    !--------------------------------------------------------------------------- 
    ! [4] Set the parameters
    !--------------------------------------------------------------------------- 

    ! initiate parameters (only touch id_pastix%iparm and id_pastix%dparm) 

    id_pastix%data = 0
    A => id_pastix%M

    Allocate(id_pastix%perm(A%n),stat=info)
    CHCKASSRT(info == 0, info)
    If (info<0) Return
    
    Allocate(id_pastix%invp(A%n),stat=info) 
    CHCKASSRT(info == 0, info)
    If (info<0) Return

    id_pastix%nrhs = 1
    nullify(id_pastix%rhs)
    id_pastix%iparm(IPARM_MODIFY_PARAMETER) = API_NO

    Call XMPH_pastix_exec(id_pastix,useRHS=.False.)
    Call XMPH_pastix_check_status(id_pastix, __LINE__, info )
    If (info < 0) Return

    ! customize parameters
    
    id_pastix%iparm(IPARM_VERBOSE)          = verbose
    id_pastix%iparm(IPARM_START_TASK)       = API_TASK_INIT
    id_pastix%iparm(IPARM_END_TASK)         = API_TASK_INIT

    select case(A%sym)
    case (SM_SYM_IsSPD)
       id_pastix%iparm(IPARM_SYM)           = API_SYM_YES
       id_pastix%iparm(IPARM_FACTORIZATION) = API_FACT_LLT
    case (SM_SYM_IsSymmetric)
       id_pastix%iparm(IPARM_SYM)           = API_SYM_YES
       id_pastix%iparm(IPARM_FACTORIZATION) = API_FACT_LDLT
    case (SM_SYM_IsGeneral )
       id_pastix%iparm(IPARM_SYM)           = API_SYM_NO
       id_pastix%iparm(IPARM_FACTORIZATION) = API_FACT_LU
    end select

    id_pastix%iparm(IPARM_MATRIX_VERIFICATION) = API_NO
    id_pastix%iparm(IPARM_RHS_MAKING)          = API_RHS_B

    !--------------------------------------------------------------------------- 
    ! [4] Initialize the XMPH_pastix_data structure.
    !--------------------------------------------------------------------------- 

    A => id_pastix%M
    Call XMPH_pastix_exec(id_pastix,useRHS=.False.)
    Call XMPH_pastix_check_status(id_pastix, __LINE__, info )
    If (info < 0) Return

    ! deactivate multi threading by default
    id_pastix%iparm(IPARM_THREAD_NBR) = 1
    id_pastix%iparm(IPARM_BINDTHRD)   = API_BIND_NO

#else
    !---------------------------------------------------------------------------
    ! [*] Libpastix unavailable
    !---------------------------------------------------------------------------
    info = ERR_LIBPASTIX_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBPASTIX_UNAVAILABLE )
    Return
#endif 

  End Subroutine XMPH_pastix_set_distributed_matrix

  !> set how threads on the pastix instance(s) should be binded to cores.
  !!
  !! When explicit binding is asked, the threads of a each pastix
  !! instance are bind to unused adjacent cores. This is necessary for optimal
  !! performance as 2 pastix instances on a same node will use the same
  !! cores by default.
  !! Here is an example : 
  !!    - If we have 8 cores (1 node; each node have 2 sockets; each socket have 4 cores)
  !!    - If we have 2 pastix instances each using 4 threads
  !!    - by default, the 2 pastix instances will use the same socket
  !!      (2 threads per core on the socket)
  !!    - by explicitly bindings the threads, we enforce 1 pastix instance per socket 
  !!      (1 pastix thread per core on each socket) 
  !!    .
  !! @see direct_solver_set_multithread

  !>
  !! @param id_pastix        the pastix instance
  !!
  !! @param thread_icntl     parameter controls (5) of multi-threadings 
  !!       - 1 binding method of threads to cores 
  !!           - -1 no multitheading
  !!           - 0 automatic binding
  !!           - 1 explicit  binding
  !!           - 2 no binding
  !!       - 2 number of nodes
  !!       - 3 number of cores per node
  !!       - 4 number of threads per process
  !!       - 5 number of processes
  !!       .
  !! @param info            status of the subroutine
  !!       - 0   success
  !!       - > 0 warning 
  !!       - < 0 error
  !!       .
  Subroutine XMPH_pastix_set_multithread (id_pastix, thread_icntl, info )

    !* Modules & co *!
    Use mph_error_mod
    Implicit none
    Include 'mpif.h'

    !* Arguments *!
    Type(XMPH_pastix_t), intent(inout) :: id_pastix
    Integer              , intent(in   ) :: thread_icntl(5)
    Integer              , intent(  out) :: info 

#if HAVE_LIBPASTIX
    !* Local variables *!

    ! Constants
    Integer, Parameter :: NO_MTHREADING = -1

    ! Scalars
    Integer :: gbrank !< MPI rank with MPI_COMM_WORLD
    Integer :: gbnp   !< MPI number of processor with MPI_COMM_WORLD

    Integer :: i
    Integer :: my_core_beg 

    !     Topology: 
    !     - n: number of nodes 
    !     - c: number of cores per node
    !     - t: number of threads per process
    !     - p: number of processes
    !     - bind : 1: use explicit binding, 0: automatic binding, 2: no binding for thread
    Integer :: topo_n, topo_c, topo_t, topo_p 
    Integer :: topo_bind

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------

    ! init local variables
    info  = 0
    topo_bind = thread_icntl(1)
    topo_n    = thread_icntl(2)
    topo_c    = thread_icntl(3)
    topo_t    = thread_icntl(4)
    topo_p    = thread_icntl(5)

    ! exit early if possible
    If ( topo_bind == NO_MTHREADING ) Return

    !---------------------------------------------------------------------------
    ! [2] Set multi-threading
    !---------------------------------------------------------------------------

    !-------------------------------------------------------------------
    ! [2.1] simplest cases
    !-------------------------------------------------------------------
    id_pastix%iparm(IPARM_THREAD_NBR) = topo_t
    if ( topo_bind == 0 ) id_pastix%iparm(IPARM_BINDTHRD)  = API_BIND_AUTO
    if ( topo_bind == 2 ) id_pastix%iparm(IPARM_BINDTHRD)  = API_BIND_NO

    !-------------------------------------------------------------------  
    ! [2.2] User specify everything case
    !-------------------------------------------------------------------

    if ( topo_bind == 1 ) then
       ! set the options
       id_pastix%iparm(IPARM_BINDTHRD)  = API_BIND_TAB

       ! obtain MPI informations
       Call MPI_Comm_rank(MPI_COMM_WORLD, gbrank, info )
       ASSRT(info == MPI_SUCCESS)

       Call MPI_Comm_size(MPI_COMM_WORLD, gbnp, info ) 
       ASSRT(info == MPI_SUCCESS)

       ! Compute the bindtab
       my_core_beg = modulo(gbrank,(topo_c / topo_t)) * topo_t;

       Allocate(id_pastix%bindtab(topo_t),STAT=info )
       CHCKASSRT(info == 0, info)
       If (info < 0) Return

       Do i=1, topo_t
          id_pastix%bindtab(i) = my_core_beg + i-1
       End Do

       ! set the bindtab 
#if defined(NEW_PASTIX_VERSION)
       Call XMPH_ARITH_pastix_fortran_bindthreads &
            (id_pastix%data, topo_t, id_pastix%bindtab)
#else
       Call XMPH_ARITH_pastix_fortran_setbindtab &
            (id_pastix%data, topo_t, id_pastix%bindtab)
#endif       
Call XMPH_pastix_check_status(id_pastix, __LINE__, info)
       If (info < 0) Return

    end if

    !---------------------------------------------------------------------------
    ! [3] Finish 
    !---------------------------------------------------------------------------
    
    ! verify
    if ( ( id_pastix%iparm(IPARM_BINDTHRD) /= API_BIND_AUTO ) .and. &
         ( id_pastix%iparm(IPARM_BINDTHRD) /= API_BIND_NO   ) .and. &
         ( id_pastix%iparm(IPARM_BINDTHRD) /= API_BIND_TAB  ) ) then
       CHCKASSRT( .False., info )
       If (info < 0) Return
    Endif

    info = 0

#else
    !---------------------------------------------------------------------------
    ! [*] Libpastix unavailable
    !---------------------------------------------------------------------------
    info = ERR_LIBPASTIX_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBPASTIX_UNAVAILABLE )
    Return
#endif 

  End subroutine XMPH_pastix_set_multithread

  ! [+] routine : XMPH_pastix_get_permutation ----------------------------------
  ! 
  !> Get the permutation done
  !! @see XMPH_sds_get_permutation()
  Subroutine XMPH_pastix_get_permutation &
       (id_pastix, perm, nperm, info)


    !* Arguments *!

    Type(XMPH_pastix_t)   , Intent(in   ) :: id_pastix
    Integer              , Intent(in   ) :: nperm
    Integer              , Intent(inout) :: perm(nperm)
    Integer              , Intent(  out) :: info

    !- End of header------------------------------------------------------------

#if HAVE_LIBPASTIX
    perm(1:nperm) = id_pastix%perm(1:nperm)
    info = 0
#else
    info = ERR_LIBPASTIX_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBPASTIX_UNAVAILABLE)
    Return
#endif

  End Subroutine XMPH_pastix_get_permutation



  ! [+] routine : XMPH_pastix_get_numstats -------------------------------------
  !
  !> get numerical statistics from at sparse direct solver.
  !! @see XMPH_sds_get_numstats
  Subroutine XMPH_pastix_get_numstats (id_pastix, piv )

    !* Arguments *!

    Type(XMPH_pastix_t), Intent(in   ) :: id_pastix
    Integer(kind=8)      , Intent(  out) :: piv

    !- End of header------------------------------------------------------------
#if HAVE_LIBPASTIX
    piv = Int(id_pastix%iparm(IPARM_STATIC_PIVOTING),8)
#else
    !---------------------------------------------------------------------------
    ! [*] Libpastix unavailable
    !---------------------------------------------------------------------------
    piv = -1_8
    Call mph_log(MSG_ERROR,LIBPASTIX_UNAVAILABLE )
    Return
#endif 

  End Subroutine XMPH_pastix_get_numstats


  ! [+] routine : XMPH_pastix_get_perfstats ------------------------------------
  !
  !> get performance statistics from at sparse direct solver.
  !! @see XMPH_sds_get_perfstats
  Subroutine XMPH_pastix_get_perfstats &
       (id_pastix, flops_estielim, flops_assemb, flops_elim )

    !* Arguments *!

    Type(XMPH_pastix_t), Intent(in   ) :: id_pastix
    Real   (kind=8)      , Intent(  out) :: flops_estielim
    Real   (kind=8)      , Intent(  out) :: flops_assemb
    Real   (kind=8)      , Intent(  out) :: flops_elim


    !* Local variables *!

    Integer(kind=8), Parameter :: UNAVAILABLE = -1

    !- End of header------------------------------------------------------------

    flops_estielim = Real(UNAVAILABLE,KIND=8)
    flops_assemb = Real(UNAVAILABLE,KIND=8)
    flops_elim = Real(UNAVAILABLE,KIND=8)

  End Subroutine XMPH_pastix_get_perfstats


  ! [+] routine : XMPH_pastix_get_memstats -------------------------------------
  ! 
  !> Get memory statistics
  !! @see XMPH_sds_get_memstats()
  Subroutine XMPH_pastix_get_memstats (id_pastix, mem )

    !* Module(s) *!

    Use XMPH_sparse_matrix_mod, Only : &
         XMPH_sm_sizeof
    Use MPH_mem_mod

    !* Arguments *!

    Type(XMPH_pastix_t), Intent(in ) :: id_pastix
    Type(MPH_mem_t          ), Intent(out) :: mem


    !* Local variables *!

    Integer(kind=8), Parameter :: UNAVAILABLE = -1
    Integer(kind=8)            :: memMtxCopy 
#if HAVE_LIBPASTIX
    Integer(kind=8)            :: pidmemused
    Integer(kind=8)            :: pidmempeak
    Integer(kind=8)            :: schurmem

    !* Explicit inferface *!
    Interface
       Subroutine sds_XMPH_pastix_getmemusage( memusage )
         Integer(kind=8), Intent(out) :: memusage
       End Subroutine sds_XMPH_pastix_getmemusage
    End Interface
#endif
    !- End of header------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------

    memMtxCopy = UNAVAILABLE
    Call MPH_mem_build(mem,UNAVAILABLE,UNAVAILABLE,UNAVAILABLE,UNAVAILABLE)

#if HAVE_LIBPASTIX
    !---------------------------------------------------------------------------
    ! [2] Set values
    !---------------------------------------------------------------------------

    If ( Associated (id_pastix%M))Then
       Call MPH_mem_add2usage(mem,XMPH_sm_sizeof(id_pastix%M))
    End If

    If ( Associated (id_pastix%schur))Then
       schurmem = id_pastix%nschur * id_pastix%ldschur * XMPH_FLOATBYTESIZE
       Call MPH_mem_add2usage(mem,schurmem)
    End If

    Call mph_sds_pastix_getmemusage(pidmemused)
    Call MPH_mem_setpidusage(mem,pidmemused)

    pidmempeak = Int(id_pastix%dparm(DPARM_MEM_MAX)     ,8)
    Call MPH_mem_setpidpeak(mem,pidmempeak)
#endif 

  End Subroutine XMPH_pastix_get_memstats


  ! [+] routine : XMPH_pastix_set_schurlist ------------------------------------
  ! 
  !> Ask to compute the dense schur during the factorization
  !! @see XMPH_sds_set_schurlist()
  !! 
  !! Warning : see README about specific calls 
  subroutine XMPH_pastix_set_schurlist  &
       (id_pastix, nschurlist, schurlist, info )

    !* Modules *!

    Use mph_error_mod
    Implicit None

    !* Arguments *!

    Type(XMPH_pastix_t), Intent(inout) :: id_pastix
    MPH_INT           , Intent(in   ) :: nschurlist
    MPH_INT           , Intent(in   ) :: schurlist(nschurlist)
    Integer              , Intent(  out) :: info 

#if HAVE_LIBPASTIX
    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] check data
    !---------------------------------------------------------------------------

    CHCKASSRT( nschurlist >= 0, info)
    If ( info < 0) Return

    !---------------------------------------------------------------------------
    ! [2] set the schur list 
    !---------------------------------------------------------------------------

    Call XMPH_ARITH_pastix_fortran_setschurunknownlist &
         (id_pastix%data, nschurlist, schurlist(1), info )
    CHCKASSRT( info >= 0, info)
    If ( info < 0) Return

    !---------------------------------------------------------------------------
    ! [3] activate the schur option 
    !---------------------------------------------------------------------------

    id_pastix%iparm(IPARM_SCHUR)               = API_YES
    id_pastix%iparm(IPARM_MATRIX_VERIFICATION) = API_NO

    ! save the leading dimension of the local schur
    id_pastix%ldschur =  nschurlist

#else
    !---------------------------------------------------------------------------
    ! [*] Libpastix unavailable
    !---------------------------------------------------------------------------
    info = ERR_LIBPASTIX_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBPASTIX_UNAVAILABLE )
    Return
#endif 

  End Subroutine XMPH_pastix_set_schurlist


  ! [+] routine : XMPH_pastix_get_schur ----------------------------------------
  ! 
  !> Get the computed dense schur during the factorization
  !! @see XMPH_sds_set_denseshur()
  !! 
  !! Warning : see doc/txt/doc_sds_pastix_mod.inc about specific calls 
  Subroutine XMPH_pastix_get_schur  &
       (id_pastix, dmschur, info )

    !* Modules *!
    Use XMPH_dense_matrix_mod
    Use mph_error_mod
#if MAPHYS_DEBUG
    Use mph_dbg_mod
#endif
    implicit none

    !* Arguments *!

    Type(XMPH_pastix_t), Intent(inout) :: id_pastix
    Type(XMPH_dense_matrix_t) , Intent(  out) :: dmschur
    Integer              , Intent(  out) :: info 

#if HAVE_LIBPASTIX
    !* Local variables *!
    Logical :: OptionSet
    Integer :: sym
    Integer :: stored

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------

    ! Check if the schur computation option was set
    OptionSet= (id_pastix%iparm(IPARM_SCHUR) == API_YES )
    CHCKASSRT( OptionSet, info  ) 
    If ( info < 0) Return

    ! Init output
    Call XMPH_dm_nullify(dmschur, info )
    CHCKASSRT( info >= 0, info )
    If ( info < 0) Return

    !---------------------------------------------------------------------------
    ! [2] Set the data
    !---------------------------------------------------------------------------

    ! Set the symmetry & storage flags
    Select Case ( id_pastix%iparm(IPARM_FACTORIZATION) )
    Case (API_FACT_LU ) ! unsymmetric
       sym    = DM_SYM_IsGeneral
       stored = DM_STORED_ALL
    Case (API_FACT_LLT ) ! SPD
       sym = DM_SYM_IsSPD
       stored = DM_STORED_LOWER
    Case (API_FACT_LDLT ) ! Symmetric
       sym = DM_SYM_IsSymmetric
       stored = DM_STORED_LOWER
    End Select

    ! Set the structure
    Call XMPH_dm_build(dmschur, &
         id_pastix%nschur, id_pastix%ldschur, id_pastix%ldschur, &
         sym,stored,id_pastix%schur, info ) 
    CHCKASSRT( info >= 0, info )
    If ( info < 0) Return

#else
    !---------------------------------------------------------------------------
    ! [*] Libpastix unavailable
    !---------------------------------------------------------------------------
    info = ERR_LIBPASTIX_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBPASTIX_UNAVAILABLE )
    Return
#endif 

  End Subroutine XMPH_pastix_get_schur


  ! [+] routine : XMPH_pastix_exec ---------------------------------------------
  !
  !> execute a pastix instance the return status of a pastix call
  !! 
  !! @param[inout] id_pastix The pastix instance to be executed
  !! @param[inout] useRHS    Use or not the pointer to the right-hand-side
  !!                         ( id_pastix%rhs ).
  Subroutine XMPH_pastix_exec(id_pastix,useRHS)

    Implicit None

    !* Argument(s) *!

    Type(XMPH_pastix_t), Intent(inout) :: id_pastix
    Logical              , Intent(in   ) :: useRHS

    !- End of header -----------------------------------------------------------

#if HAVE_LIBPASTIX
    If ( useRHS ) Then
       Call XMPH_ARITH_pastix_fortran                 &
            (id_pastix%data, id_pastix%comm, &
            id_pastix%M%n,id_pastix%M%csc,   &
            id_pastix%M%i,id_pastix%M%v,     &
            id_pastix%perm,id_pastix%invp,   &
            id_pastix%rhs ,id_pastix%nrhs,   & 
            id_pastix%iparm,id_pastix%dparm )
    Else
       Call XMPH_ARITH_pastix_fortran                 &
            (id_pastix%data, id_pastix%comm, &
            id_pastix%M%n,id_pastix%M%csc,   &
            id_pastix%M%i,id_pastix%M%v,     &
            id_pastix%perm,id_pastix%invp,   &
            id_pastix%dummy_rhs,id_pastix%nrhs, & 
            id_pastix%iparm,id_pastix%dparm )
    End If
#endif

  End Subroutine XMPH_pastix_exec


  ! [+] routine : XMPH_pastix_check_status -------------------------------------
  !
  !> check the return status of a pastix call
  !! 
  !! @param[in ] id_pastix The pastix instance to be checked
  !! @param[in ] line      The line where the check was performed 
  !! @param[out] info      The maphys return status (0< error, 0: ok, 1> warn)
  !!
  Subroutine XMPH_pastix_check_status(id_pastix, line, info )

    !* Module(s) *!
    Use mph_log_mod
    Implicit None

    !* Argument(s) *!

    Type(XMPH_pastix_t), Intent(in   ) :: id_pastix
    Integer              , Intent(in   ) :: line
    Integer              , Intent(  out) :: info

#if HAVE_LIBPASTIX
    !* Local variables *!
    
    ! strings
    Character(len=MAPHYS_STRL) :: msg
    
    !- End of header -----------------------------------------------------------

    ! Get return code 
    info = id_pastix%iparm(IPARM_ERROR_NUMBER);
    
    ! On success, nothing todo
    If (info == 0) Return

    ! Print a maximum of information, to help debuging.
    Write(msg,'(2A,I6,A,I10)') &
         Trim(FLNAME),":line",line, &
         ". Details : iparm(IPARM_ERROR_NUMBER)=", info
    Call mph_log(MSG_ERROR,Trim(msg))

    ! pastix only have errors
    If ( info > 0 ) info = -info

#else
    !---------------------------------------------------------------------------
    ! [*] Libpastix unavailable
    !---------------------------------------------------------------------------
    info = ERR_LIBPASTIX_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBPASTIX_UNAVAILABLE )
    Return
#endif 

  End Subroutine XMPH_pastix_check_status

  ! [+] routine : XMPH_pastix_sparse_matrix_check ----------------------------
  !
  !> maphys driver for pastix checkmatrix.
  !!
  !! The driver symmetrize the matrix structure (in CSC)
  !!
  !! @param[in,out] M      The matrix to be checked
  !! @param[in    ] comm   The mpi communicator
  !! @param[   out] info   The return status
  subroutine XMPH_pastix_sparse_matrix_check(M, comm, info )
    
    !* Module(s) *!
    Use XMPH_sparse_matrix_mod
    Use mph_error_mod
#if MAPHYS_DEBUG
    Use mph_dbg_mod
#endif
    implicit none

    !* Argument(s) *!
    Type(XMPH_sparse_matrix_t), intent(inout) :: M        
    Integer              , intent(inout) :: comm
    Integer              , intent(  out) :: info 

#if HAVE_LIBPASTIX
    !* Local variables *!

    ! constants
#if MAPHYS_XMPH_PASTIX_VERBOSE
    Integer, Parameter                   :: verbose = API_VERBOSE_YES
#else
    Integer, Parameter                   :: verbose = API_VERBOSE_NOT
#endif

    ! scalars
    Integer                              :: XMPH_pastix_sym
    Integer                              :: allow_size_modif

    ! type
    pastix_data_ptr_t                    :: check_data
    
    !- End of header -----------------------------------------------------------
    
    !---------------------------------------------------------------------------
    ! [1] Preprocess the matrix in the right format
    !---------------------------------------------------------------------------

    ! On symmetric matrices, 
    ! only the inferior triangular part must be given to Pastix

    Select Case (M%sym)
    Case( SM_SYM_IsSPD, SM_SYM_IsSymmetric) 
       Call XMPH_sm_transposeUpper(M)
    Case (SM_SYM_IsGeneral)
       Continue
    Case Default
       Continue
    End Select

    ! Convert to csc if neccessary

    Call XMPH_sm_convert(M, SM_FMT_CSC, info )
    CHCKASSRT(info >= 0, info)
    If (info < 0) Return

    !---------------------------------------------------------------------------
    ! [2] Select strategy according to the symetry
    !---------------------------------------------------------------------------

    If (M%sym == SM_SYM_IsSPD      ) XMPH_pastix_sym = API_SYM_YES
    If (M%sym == SM_SYM_IsSymmetric) XMPH_pastix_sym = API_SYM_YES
    If (M%sym == SM_SYM_IsGeneral  ) XMPH_pastix_sym = API_SYM_NO

    if (M%sym == SM_SYM_IsSPD      ) allow_size_modif = API_NO
    If (M%sym == SM_SYM_IsSymmetric) allow_size_modif = API_NO
    If (M%sym == SM_SYM_IsGeneral  ) allow_size_modif = API_YES

    !---------------------------------------------------------------------------
    ! [3] check the matrix 
    !---------------------------------------------------------------------------

    Call XMPH_ARITH_pastix_fortran_checkmatrix &
         &(check_data, comm, verbose,     &
         & XMPH_pastix_sym, allow_size_modif,  & 
         & M%n, M%csc, M%i, M%v, &
         & -1, 1)

    !---------------------------------------------------------------------------
    ! [4] If matrix structure need to be symmetrized, reallocate & call pastix.
    !---------------------------------------------------------------------------

    !  resize the matrix (add symetric elements)
    if (M%csc(M%n+1) - 1 /= M%nnz ) then

       Deallocate(M%i)
       Deallocate(M%j)
       Deallocate(M%v)

       M%nnz = M%csc(M%n+1)-1
       Allocate(M%i(M%nnz),STAT=info )
       CHCKASSRT(info == 0, info)
       If (info < 0) Return

       Allocate(M%j(M%nnz),STAT=info )
       CHCKASSRT(info == 0, info)
       If (info < 0) Return

       Allocate(M%v(M%nnz),STAT=info )
       CHCKASSRT(info == 0, info)
       If (info < 0) Return

       Call XMPH_ARITH_pastix_fortran_checkmatrix_end &
            &(check_data, verbose, M%i, M%v, 1)
    Endif

    ! rebuild the j array of sparse_matrix (maintain a valid i,j,v)

    Call XMPH_sm_ptr2ind(M%n+1,M%csc,M%j)

#else
    !---------------------------------------------------------------------------
    ! [*] Libpastix unavailable
    !---------------------------------------------------------------------------
    info = ERR_LIBPASTIX_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBPASTIX_UNAVAILABLE )
    Return
#endif 

  End Subroutine XMPH_pastix_sparse_matrix_check


  ! [+] routine : XMPH_pastix_free_schur -----------------------------------------
  ! 
  !> Free the schur complement
  !! @see XMPH_sds_free_schur()
  Subroutine XMPH_pastix_free_schur (id_pastix, info )
    Implicit None

    !* Arguments *!

    Type(XMPH_pastix_t)        , Intent(inout) :: id_pastix
    Integer              , Intent(  out) :: info

    !- End of header------------------------------------------------------------

#if HAVE_LIBPASTIX
    info = MPH_SUCCESS
    If (.Not. Associated(id_pastix%SCHUR))Then
       info = 1
       Return
    End If
    
    Deallocate(id_pastix%SCHUR)
#else
    !---------------------------------------------------------------------------
    ! [*] Libpastix unavailable
    !---------------------------------------------------------------------------
    info = ERR_LIBPASTIX_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBPASTIX_UNAVAILABLE )
    Return
#endif 


  End Subroutine XMPH_pastix_free_schur



End Module XMPH_sds_pastix_mod
