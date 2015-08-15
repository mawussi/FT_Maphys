! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"

! define the partitioner to use (METIS by default)
#ifndef MAPHYS_MUMPS_ICNTL7
#define MAPHYS_MUMPS_ICNTL7 5
#endif

! define if MUMPS should log its messages (default : No)
#ifndef MAPHYS_MUMPS_LOG
#define MAPHYS_MUMPS_LOG 0
#endif

! [+] module : XMPH_sds_mumps_mod ----------------------------------------------
!
!> MaPHyS MUMPS adapter module
!! @see XMPH_sds_mod
!! @see doc_sds_mumps_mod.inc
Module XMPH_sds_mumps_mod

  !* Module(s) *!
  Use XMPH_sds_mumps_type, Only : &
       XMPH_mumps_t  => XMPH_ARITHmumps_struc ! renamed structure.

#if HAVE_LIBMUMPS

  !* No implicit typing *!
  Implicit None

  !* Private constants *!
  Character(len=MAPHYS_STRL), Parameter, Private :: &
       FLNAME = "XMPH_ARITHmph_sds_mumps_mod.F90"

  !* Explicit interfaces *! 
  Interface x_mumps
     !**** 
     Subroutine dmumps(id)
       Include 'dmumps_struc.h'
       Type(DMUMPS_STRUC) :: id
     End Subroutine dmumps
     !****
     Subroutine smumps(id)
       Include 'smumps_struc.h'
       Type(SMUMPS_STRUC) :: id
     End Subroutine smumps
     !****
     Subroutine cmumps(id)
       Include 'cmumps_struc.h'
       type(CMUMPS_STRUC) :: id
     End subroutine cmumps
     !****
     Subroutine zmumps(id)
       include 'zmumps_struc.h'
       type(ZMUMPS_STRUC) :: id
     End subroutine zmumps
     !****
  End Interface

#else
  Use mph_log_mod
  Implicit None
  Integer, Parameter, Private :: ERR_LIBMUMPS_UNAVAILABLE = -9999
  Character(len=MAPHYS_STRL), Parameter, Private :: &
       LIBMUMPS_UNAVAILABLE = "library libXMPH_ARITHmumps is unavailable"
#endif


  !* Access specifiers *!
  Public :: XMPH_mumps_available

  Private :: XMPH_mumps_check_status
  Private :: XMPH_mumps_nullify_user_arrays 

  Public :: XMPH_mumps_set_MPIcommunicator 
  Public :: XMPH_mumps_set_matrix
  Public :: XMPH_mumps_set_distributed_matrix 
  Public :: XMPH_mumps_analyze 
  Public :: XMPH_mumps_set_denseschur
  Public :: XMPH_mumps_factorize 
  Public :: XMPH_mumps_solve_RHS
  Public :: XMPH_mumps_solve_sparseRHS
  Public :: XMPH_mumps_finalize 
  Public :: XMPH_mumps_get_permutation
  Public :: XMPH_mumps_get_numstats
  Public :: XMPH_mumps_get_memstats
  Public :: XMPH_mumps_get_perfstats
  Public :: XMPH_mumps_set_schurlist
  Public :: XMPH_mumps_get_schur
  Public :: XMPH_mumps_free_schur

  !* Routines *!

Contains

  ! [+] function : XMPH_mumps_available ----------------------------------------
  ! 
  !> Return 0 if mumps library is available or -1 if not.
  Function XMPH_mumps_available()
    Integer :: XMPH_mumps_available
#if HAVE_LIBMUMPS    
    XMPH_mumps_available = 0
#else
    XMPH_mumps_available = -1
#endif
  End Function XMPH_mumps_available

  ! [+] routine : XMPH_mumps_check_status --------------------------------------
  ! 
  !> check the return status of the mumps structure
  !! 
  !! @param[in ] id_mumps The mumps instance to be checked
  !! @param[in ] line     The line where the check was performed 
  !! @param[out] info     The maphys return status (0< error, 0: ok, 1> warn)
  Subroutine XMPH_mumps_check_status(id_mumps, line, info )

    !* Module(s) *!
    Use mph_log_mod
    Implicit None

    !* Argument(s) *!

    Type(XMPH_mumps_t)        , Intent(in   ) :: id_mumps
    Integer              , Intent(in   ) :: line
    Integer              , Intent(  out) :: info

#if HAVE_LIBMUMPS    

    !* Local variables *!

    ! scalars
    Integer :: msgclass
    
    ! strings
    Character(len=MAPHYS_STRL) :: msg
    
    !- End of header -----------------------------------------------------------

    ! Set return code
    info =  id_mumps%infog(1)

    ! On success, nothing todo
    If (id_mumps%infog(1) == 0) Return

    ! Set if it is a warning or an error
    If (id_mumps%infog(1) >  0) msgclass = MSG_WARNING
    If (id_mumps%infog(1) <  0) msgclass = MSG_ERROR

    ! Print a maximum of information, to help debuging.
    Write(msg,'(2A,I6,A,I10,A,I10,A,I10,A,I10)') &
         Trim(FLNAME),":line",line, &
         ". Details : infog(1)=", id_mumps%infog(1),& 
         "; infog(2)=", id_mumps%infog(2),& 
         "; info (1)=", id_mumps%info(1),& 
         "; info (2)=", id_mumps%info(2)

    Call mph_log(msgclass,Trim(msg))

#else
    info = ERR_LIBMUMPS_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBMUMPS_UNAVAILABLE)
    Return
#endif 

  End Subroutine XMPH_mumps_check_status

  ! [+] routine : XMPH_mumps_nullify_user_arrays -------------------------------
  !
  !> Nullify the arrays set by the user in the mumps structure.
  !! 
  Subroutine XMPH_mumps_nullify_user_arrays (id_mumps)

    !* Argument(s) *!
    
    Type(XMPH_mumps_t), intent(inout) :: id_mumps

    !- End of header -----------------------------------------------------------

#if HAVE_LIBMUMPS

    Nullify(id_mumps%IRN)
    Nullify(id_mumps%JCN)
    Nullify(id_mumps%A)
    Nullify(id_mumps%IRN_loc)
    Nullify(id_mumps%JCN_loc)
    Nullify(id_mumps%A_loc  )
    Nullify(id_mumps%LISTVAR_SCHUR)
    Nullify(id_mumps%SCHUR)
    Nullify(id_mumps%RHS)
    Nullify(id_mumps%RHS_SPARSE)
    Nullify(id_mumps%IRHS_SPARSE)
    Nullify(id_mumps%IRHS_PTR)

#endif 

  End Subroutine XMPH_mumps_nullify_user_arrays

  ! [+] routine : XMPH_mumps_set_MPIcommunicator -------------------------------
  ! 
  !> set the MPI communicator to use.
  !! @see XMPH_sds_set_MPIcommunicator()
  Subroutine XMPH_mumps_set_MPIcommunicator (id_mumps, comm, info )

    !* Module(s) *!

    Use mph_log_mod
    Use mph_error_mod
    Implicit none
    Include 'mpif.h'

    !* Argument(s) *!

    Type(XMPH_mumps_t), intent(inout) :: id_mumps
    Integer      , intent(in   ) :: comm
    Integer      , intent(  out) :: info

#if HAVE_LIBMUMPS

    !- End of header -----------------------------------------------------------

    Call MPI_Comm_dup(comm, id_mumps%COMM , info )
    ASSRT( info == MPI_SUCCESS )

    Call XMPH_mumps_nullify_user_arrays(id_mumps)
    info = 0

#else
    info = ERR_LIBMUMPS_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBMUMPS_UNAVAILABLE)
    Return
#endif 



  End Subroutine XMPH_mumps_set_MPIcommunicator

  ! [+] routine : XMPH_mumps_set_matrix ----------------------------------------
  !
  !> Set the matrix centralized on process 0.
  !! @see XMPH_sds_set_matrix()
  subroutine XMPH_mumps_set_matrix &
       (id_mumps, M, info )

    !* Module(s) *!

    Use XMPH_sparse_matrix_mod
    Use mph_error_mod
    Use mph_log_mod
    Implicit None
    Include "mpif.h"

    !* Argument(s) *!

    type(XMPH_mumps_t)        , intent(inout) :: id_mumps
    type(XMPH_sparse_matrix_t), intent(in   ) :: M
    integer              , intent(  out) :: info

#if HAVE_LIBMUMPS

    !* Local variables *!

    ! scalars
    Integer :: rank
    Logical :: maySegfault
#if MAPHYS_MUMPS_LOG
    Integer :: verbosity
#endif

    !- End of header------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1]  Get default parameters
    !---------------------------------------------------------------------------

    ! perform some checking which may lead to segfault.
    ! If pastix runs on near-dense matrix and the
    ! number of elements is huge, segfault may occurs.
    maySegfault =( M%m*M%n /= M%nnz ).And.(M%nnz >= Huge(M%nnz)/(MPH_INTKIND))
    If( maySegfault )&
         Call MPH_LogWithInfo(MSG_WARNING,M%nnz,&
         "(mumps wrapper) matrix is not dense,&
         & and its nnz is near the limit of the integer.&
         & Segmentation fault may occur. nnz =")

    ! 
    id_mumps%job = -1
    id_mumps%par = 1

    Select Case (M%sym)
    Case (SM_SYM_IsGeneral  );  id_mumps%sym = 0
    Case (SM_SYM_IsSPD      );  id_mumps%sym = 1
    Case (SM_SYM_IsSymmetric);  id_mumps%sym = 2
    End Select

    Call x_mumps(id_mumps)
    Call XMPH_mumps_check_status(id_mumps,__LINE__,info)
    If (info < 0) Return

    !---------------------------------------------------------------------------
    ! [2]  Set the parameters (icntl/cntl)
    !---------------------------------------------------------------------------

#if MAPHYS_MUMPS_LOG
    Call mph_log_GetUnit(MSG_ERROR  , id_mumps%ICNTL(1))
    Call mph_log_GetUnit(MSG_WARNING, id_mumps%ICNTL(2))
    Call mph_log_GetUnit(MSG_STD    , id_mumps%ICNTL(3))
    Call mph_log_GetVerbosity(verbosity)
    If (verbosity >= MSG_ERROR   ) id_mumps%ICNTL(4) = 1
    If (verbosity >= MSG_DEBUG   ) id_mumps%ICNTL(4) = 2
#else
    id_mumps%ICNTL(1) = 0
    id_mumps%ICNTL(2) = 0
    id_mumps%ICNTL(3) = 0
    id_mumps%ICNTL(4) = 0
#endif

    ! others
    id_mumps%ICNTL(5) = 0   ! the matrix is given is assembled format
    id_mumps%ICNTL(18)= 0   ! the matrix is centralised 

    ! Pivoting is activated by default
    ! TODO : need to set these properties through an option in ICNTL...
    id_mumps%ICNTL(14)=  20  ! set the estimated memory increase 
    id_mumps%CNTL(1)  = 1.d-9 ! deactivate pivoting

    !---------------------------------------------------------------------------
    ! [3]  Set the matrix
    !---------------------------------------------------------------------------

    Call MPI_Comm_rank(id_mumps%comm, rank, info)
    ASSRT( info == MPI_SUCCESS )
    
    If(rank == 0)Then
       id_mumps%N    = M%m
       id_mumps%NZ   = M%nnz
       id_mumps%IRN  => M%i
       id_mumps%JCN  => M%j
       id_mumps%A    => M%v
    Endif

    ! set return status
    If (maySegfault) info = 1

#else
    info = ERR_LIBMUMPS_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBMUMPS_UNAVAILABLE)
    Return
#endif 


  End Subroutine XMPH_mumps_set_matrix

  ! [+] routine : XMPH_mumps_set_distributed_matrix ----------------------------
  ! 
  !> set the matrix which is distributed across mumps%comm
  !! @see XMPH_sds_set_distributed_matrix()
  Subroutine XMPH_mumps_set_distributed_matrix (id_mumps, M, global_n, info )

    !* Module(s) *!

    Use XMPH_sparse_matrix_mod
    Implicit None

    !* Argument(s) *!

    Type(XMPH_mumps_t)        , Intent(inout) :: id_mumps        
    Type(XMPH_sparse_matrix_t), Intent(in   ) :: M       
    Integer              , Intent(in   ) :: global_n 
    Integer              , Intent(  out) :: info    

#if HAVE_LIBMUMPS

    !* Local variable(s) *!

    ! constants
    Integer, Parameter :: verbosity = 0 ! 0 (no), 6 (max)

    !- End of header----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Initialize the instance
    !---------------------------------------------------------------------------

    id_mumps%JOB = -1
    id_mumps%SYM =  M%sym
    id_mumps%PAR =  1
    Call x_mumps(id_mumps)
    Call XMPH_mumps_check_status(id_mumps,__LINE__,info)
    If (info < 0) Return

    id_mumps%ICNTL(3) = verbosity

    !---------------------------------------------------------------------------
    ! [2] Give the distributed matrix 
    !---------------------------------------------------------------------------

    id_mumps%ICNTL(18) = 3
    id_mumps%ICNTL(14) = 70

    id_mumps%N       =  global_n
    id_mumps%NZ_loc  =  M%nnz
    id_mumps%IRN_loc => M%i
    id_mumps%JCN_loc => M%j
    id_mumps%A_loc   => M%v

#else
    info = ERR_LIBMUMPS_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBMUMPS_UNAVAILABLE)
    Return
#endif 


  End Subroutine XMPH_mumps_set_distributed_matrix

  ! [+] routine : XMPH_mumps_analyze -------------------------------------------
  ! 
  !> Analyze the linear system.
  !! @see XMPH_sds_analyze()
  Subroutine XMPH_mumps_analyze (id_mumps, info )

    Implicit None

    !* Argument(s) *!

    Type(XMPH_mumps_t), Intent(inout) :: id_mumps
    Integer      , Intent(  out) :: info

#if HAVE_LIBMUMPS

    !- End of header -----------------------------------------------------------

    id_mumps%ICNTL(7) = MAPHYS_MUMPS_ICNTL7
    id_mumps%job = 1 
    Call x_mumps(id_mumps)
    Call XMPH_mumps_check_status(id_mumps,__LINE__,info)
    If (info < 0) Return

#else
    info = ERR_LIBMUMPS_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBMUMPS_UNAVAILABLE)
    Return
#endif 


  End Subroutine XMPH_mumps_analyze

  ! [+] routine : XMPH_mumps_set_denseschur ------------------------------------
  ! 
  !> Ask to compute the dense schur
  !! @see XMPH_sds_set_denseshur()
  Subroutine XMPH_mumps_set_denseschur  &
       (id_mumps, schurlist, nschurlist, schur, nschur, ldschur, info)

    !* Modules *!
    Use mph_error_mod
    Implicit None

    !* Arguments *!

    Type(XMPH_mumps_t)           , Intent(inout) :: id_mumps
    Integer                 , Intent(in   ) :: nschurlist
    Integer                 , Intent(in   ) :: schurlist(nschurlist)
    Integer                 , Intent(in   ) :: nschur,ldschur

    XMPH_FLOAT  , Pointer , Intent(inout) :: schur(:)
    Integer                 , Intent(  out) :: info

#if HAVE_LIBMUMPS

    !* Local variables
    MPH_INT :: i

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Check arguments, initialize variables, etc.
    !---------------------------------------------------------------------------

    ! check

    Call mph_check( nschurlist > 0, __LINE__, info, FLNAME)
    If ( info < 0) Return

    Call mph_check( nschur     > 0, __LINE__, info, FLNAME)
    If ( info < 0) Return

    Call mph_check( ldschur    > 0, __LINE__, info, FLNAME)
    If ( info < 0) Return

    Call mph_check( Associated (schur), __LINE__, info, FLNAME)
    If ( info < 0) Return

    ! alloc
    Allocate( id_mumps%LISTVAR_SCHUR(nschurlist), STAT=info )
    Call mph_check( info == 0 , __LINE__, info, FLNAME )
    If ( info < 0) Return

    !---------------------------------------------------------------------------
    ! [2] Set mumps arguments
    !---------------------------------------------------------------------------
    ! warning argument "ldschur" is ignored here

    ! activate the option 
    id_mumps%ICNTL(19) = 1

    ! set input (list of variable)
    Do i=1, nschurlist
       id_mumps%LISTVAR_SCHUR(i) = schurlist(i)
    End Do

    ! set output (schur)
    id_mumps%SIZE_SCHUR = nschur
    id_mumps%SCHUR => schur

#else
    info = ERR_LIBMUMPS_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBMUMPS_UNAVAILABLE)
    Return
#endif 

  End Subroutine XMPH_mumps_set_denseschur

  ! [+] routine : XMPH_mumps_factorize -----------------------------------------
  ! 
  !> Factorize the matrix 
  !! @see XMPH_sds_factorize()
  Subroutine XMPH_mumps_factorize (id_mumps, info )

    !* Module(s) *!
    Use mph_error_mod
    Use XMPH_dense_matrix_mod
    Implicit None

    !* Argument(s) *!

    Type(XMPH_mumps_t), intent(inout) :: id_mumps
    Integer      , intent(  out) :: info

    !* Local variable(s) *!

    ! scalars
    Logical :: DistSchurAsked
    Logical :: CentSchurAsked
    MPH_INT :: i,j, nschur

    ! arrays
    XMPH_FLOAT, pointer :: schur (:)

    !- End of header -----------------------------------------------------------

#if HAVE_LIBMUMPS    

    !---------------------------------------------------------------------------
    ! [1] Init 
    !---------------------------------------------------------------------------

    ! Set a few flags
    CentSchurAsked = ( id_mumps%ICNTL(19) == 1 )
    DistSchurAsked = ( id_mumps%ICNTL(19) == 2 ).OR.( id_mumps%ICNTL(19) == 3 )

    ! Allocate the distributed schur complement
    If ( DistSchurAsked )Then
       id_mumps%SCHUR_LLD = id_mumps%SCHUR_MLOC
       i = id_mumps%SCHUR_LLD * ( id_mumps%SCHUR_NLOC - 1) + &
            id_mumps%SCHUR_MLOC
       Allocate( id_mumps%SCHUR(i), STAT=info )
       CHCKASSRT( info == 0, info )
       If ( info < 0) Return
    End If

    !---------------------------------------------------------------------------
    ! [1] Perform factorization 
    !---------------------------------------------------------------------------

    id_mumps%job = 2 
    Call x_mumps(id_mumps)
    Call XMPH_mumps_check_status(id_mumps,__LINE__,info)
    If (info < 0) Return

    !---------------------------------------------------------------------------
    ! [2] Post-treat the centralized schur complement if it was asked
    !---------------------------------------------------------------------------
    ! we want the complete schur stored by columns ( not by rows )

    If ( CentSchurAsked ) Then
       schur => id_mumps%schur
       nschur = id_mumps%SIZE_SCHUR
       unsymmetric_matrix:If (id_mumps%sym == 0) Then
          Call XMPH_transpose_matrix_1D (schur,nschur)
       Else 
          ! complete the upper triangular part
          Do i=1,nschur
             Do j=1,i-1
                schur((j-1)*nschur+i) = schur((i-1)*nschur+j)
             End Do
          End Do
       Endif unsymmetric_matrix

       Nullify(schur)
    Endif 

#else
    info = ERR_LIBMUMPS_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBMUMPS_UNAVAILABLE)
    Return
#endif 


  End Subroutine XMPH_mumps_factorize

  ! [+] routine : XMPH_mumps_solve_RHS -----------------------------------------
  ! 
  !> Solve the linear system with a dense, centralized right-hand-side
  !! @see XMPH_sds_solve_RHS()
  Subroutine XMPH_mumps_solve_RHS  &
       (id_mumps, rhs, nrhs, ldrhs, info )

    !* Module(s) *!

    Use mph_log_mod
    Use mph_error_mod
    Implicit None
    Include "mpif.h"

    !* Argument(s) *!

    Type(XMPH_mumps_t)          , intent(inout) :: id_mumps
    Integer                , intent(in   ) :: nrhs,ldrhs
    XMPH_FLOAT, target   , intent(inout) :: rhs(ldrhs*nrhs)
    Integer                , intent(  out) :: info

    !* Local variable(s) *!

    ! scalars
    Integer               :: rank

    !- End of header -----------------------------------------------------------

#if HAVE_LIBMUMPS    

    !---------------------------------------------------------------------------
    ! [1] Set the controls
    !---------------------------------------------------------------------------

    ! solve
    id_mumps%job = 3 
    ! on a centralised dense right and side
    id_mumps%ICNTL(20) =  0 
    ! on a centralised dense solution
    id_mumps%ICNTL(21) =  0

    !---------------------------------------------------------------------------
    ! [2] Set the RHS
    !---------------------------------------------------------------------------

    ! given by root = 0
    Call MPI_Comm_rank(id_mumps%comm, rank, info)
    Call mph_assert(info == MPI_SUCCESS, __LINE__, FLNAME )
    if ( rank == 0 ) then 
       id_mumps%rhs  => rhs
       id_mumps%nrhs =  nrhs
       id_mumps%lrhs =  ldrhs
    end if

    !---------------------------------------------------------------------------
    ! [3] Perform solve step
    !---------------------------------------------------------------------------

    Call x_mumps(id_mumps)
    Call XMPH_mumps_check_status(id_mumps,__LINE__,info)
#endif 

  End Subroutine XMPH_mumps_solve_RHS

  ! [+] routine : XMPH_mumps_solve_sparseRHS -----------------------------------
  ! 
  !> Solve the linear system with a sparse right-hand-side in coordinate format
  !! @see XMPH_sds_solve_sparseRHS()
  Subroutine XMPH_mumps_solve_sparseRHS &
       (id_mumps, nx, ldx, x, sp_rhs, info)

    !* Module(s) *!
    Use XMPH_sparse_matrix_mod
    Use mph_error_mod
    Implicit None

    !* Argument(s) *!

    type(XMPH_mumps_t), intent(inout) :: id_mumps
    Integer      , intent(in   ) :: nx,ldx
    XMPH_FLOAT , target , intent(inout) :: x(nx*ldx)

    ! warning in CSC format
    Type(XMPH_sparse_matrix_t), intent(in   ) :: sp_rhs
    Integer              , intent(  out) :: info

    !* Local variable(s) *!

    Integer, Parameter :: verbosity = 0 ! 0 (no), 6 (max)

    !- End of header -----------------------------------------------------------

#if HAVE_LIBMUMPS    

    !---------------------------------------------------------------------------
    ! [1] Check arguments
    !---------------------------------------------------------------------------

    Call mph_check(ldx == sp_rhs%m, __LINE__, info, FLNAME)
    If ( info < 0) Return
    Call mph_check(nx  == sp_rhs%n, __LINE__, info, FLNAME)
    If ( info < 0) Return

    !---------------------------------------------------------------------------
    ! [2] Set controls
    !---------------------------------------------------------------------------

    id_mumps%ICNTL(3)  =  verbosity

    ! solve
    id_mumps%JOB       =  3
    ! a sparse RHS option
    id_mumps%ICNTL(20) =  1 

    !---------------------------------------------------------------------------
    ! [3] Set the RHS/SOL
    !---------------------------------------------------------------------------

    ! the sparse RHS
    id_mumps%NRHS      =  sp_rhs%n
    id_mumps%NZ_RHS    =  sp_rhs%nnz

    id_mumps%IRHS_SPARSE => sp_rhs%i
    id_mumps%RHS_SPARSE  => sp_rhs%v
    id_mumps%IRHS_PTR    => sp_rhs%csc

    ! the solution
    id_mumps%LRHS      =  ldx
    id_mumps%RHS       => x

    !---------------------------------------------------------------------------
    ! [4] Perform solve step
    !---------------------------------------------------------------------------

    Call x_mumps(id_mumps)
    Call XMPH_mumps_check_status(id_mumps,__LINE__,info)

#else
    info = ERR_LIBMUMPS_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBMUMPS_UNAVAILABLE)
    Return
#endif 

  End Subroutine XMPH_mumps_solve_sparseRHS

  ! [+] routine : XMPH_mumps_finalize ------------------------------------------
  ! 
  !> Free the mumps instance
  !! @see XMPH_sds_finalize()
  Subroutine XMPH_mumps_finalize (id_mumps, info )
    
    Implicit None

    !* Argument(s) *!

    Type(XMPH_mumps_t), Intent(inout) :: id_mumps
    Integer      , Intent(  out) :: info

    !- End of header -----------------------------------------------------------

#if HAVE_LIBMUMPS    

    ! Deallocate arrays allocated by the wrapper for the schur 
    Select Case( id_mumps%ICNTL(19) )
    Case(1)      ! centralized schur
       Deallocate(id_mumps%LISTVAR_SCHUR)
    Case(2,3)    ! distributed schur
       Deallocate(id_mumps%LISTVAR_SCHUR)
       Deallocate(id_mumps%SCHUR)
    Case Default ! no schur
       Continue
    End Select

    ! Deassociate arrays set by the wrapper
    Call XMPH_mumps_nullify_user_arrays(id_mumps)

    ! Destroy the instance
    id_mumps%job=-2
    Call x_mumps(id_mumps)
    Call XMPH_mumps_check_status(id_mumps,__LINE__,info)

#else
    info = ERR_LIBMUMPS_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBMUMPS_UNAVAILABLE)
    Return
#endif 

  End Subroutine XMPH_mumps_finalize

  ! [+] routine : XMPH_mumps_get_permutation -----------------------------------
  ! 
  !> Get the permutation done
  !! @see XMPH_sds_get_permutation()
  Subroutine XMPH_mumps_get_permutation &
       (id_mumps, perm, nperm, info)


    !* Arguments *!

    Type(XMPH_mumps_t)   , Intent(in   ) :: id_mumps
    Integer              , Intent(in   ) :: nperm
    Integer              , Intent(inout) :: perm(nperm)
    Integer              , Intent(  out) :: info

    !- End of header------------------------------------------------------------

#if HAVE_LIBMUMPS
    perm(1:nperm) = id_mumps%SYM_PERM(1:nperm)
    info = 0
#else
    info = ERR_LIBMUMPS_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBMUMPS_UNAVAILABLE)
    Return
#endif

  End Subroutine XMPH_mumps_get_permutation


  ! [+] routine : XMPH_mumps_get_numstats --------------------------------------
  ! 
  !> Get numerical statistics
  !! @see XMPH_sds_get_numstats()
  Subroutine XMPH_mumps_get_numstats (id_mumps, piv )

    !* Arguments *!

    Type(XMPH_mumps_t)        , Intent(in   ) :: id_mumps
    Integer(kind=8)      , Intent(  out) :: piv

    !* Local variables *!

    Integer(kind=8), Parameter :: UNAVAILABLE = -1

    !- End of header------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------

    piv = UNAVAILABLE

    !---------------------------------------------------------------------------
    ! [2] Get number of static pivots
    !---------------------------------------------------------------------------

#if HAVE_LIBMUMPS    
    If ( id_mumps%JOB >= 2 ) piv = id_mumps%infog(25)
#endif 

  End Subroutine XMPH_mumps_get_numstats

  ! [+] routine : XMPH_mumps_get_memstats --------------------------------------
  ! 
  !> Get several statistics
  !! @see XMPH_sds_get_memstats()
  Subroutine XMPH_mumps_get_memstats (id_mumps, mem )

    !* Module(s) *!

    Use MPH_mem_mod

    !* Arguments *!

    Type(XMPH_mumps_t)        , Intent(in   ) :: id_mumps
    Type(MPH_mem_t          ), Intent(out) :: mem

    !* Local variables *!

    Integer(kind=8), Parameter :: UNAVAILABLE = -1
#if HAVE_LIBMUMPS    
    Integer(kind=8) :: memused
    Integer(kind=8) :: mempeak
#endif

    !- End of header------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------

    Call MPH_mem_build(mem,UNAVAILABLE,UNAVAILABLE,UNAVAILABLE,UNAVAILABLE)

#if HAVE_LIBMUMPS    

    !---------------------------------------------------------------------------
    ! [2] Set the memory used
    !---------------------------------------------------------------------------

    If ( id_mumps%JOB >= 2 )Then

       ! count internal memory of id_mumps
       mempeak = 1.0E6 * id_mumps%infog(19) 
       memused = 1.0E6 * id_mumps%infog(22) 

       Call MPH_mem_setusage(mem,memused)
       Call MPH_mem_setpeak (mem,mempeak)
       Call MPH_mem_setpidpeak  (mem,0_8)
       Call MPH_mem_setpidusage (mem,0_8)

       ! count the allocated schur
       If ((id_mumps%ICNTL(19) == 2).OR.(id_mumps%ICNTL(19) == 3))Then
          memused = id_mumps%SCHUR_LLD * ( id_mumps%SCHUR_NLOC - 1) + &
            id_mumps%SCHUR_MLOC
          memused = memused * XMPH_FLOATBYTESIZE
          Call MPH_mem_add2usage(mem,memused)
       End If

    End If

#endif 

  End Subroutine XMPH_mumps_get_memstats

  ! [+] routine : XMPH_mumps_get_perfstats -------------------------------------
  ! 
  !> Get performance statistics
  !! @see XMPH_sds_get_perfstats()
  Subroutine XMPH_mumps_get_perfstats &
       (id_mumps, flops_estielim, flops_assemb, flops_elim )
    !* Arguments *!

    Type(XMPH_mumps_t)        , Intent(in   ) :: id_mumps
    Real   (kind=8)      , Intent(  out) :: flops_estielim
    Real   (kind=8)      , Intent(  out) :: flops_assemb
    Real   (kind=8)      , Intent(  out) :: flops_elim

    !* Local variables *!
    Integer(kind=8), Parameter :: UNAVAILABLE = -1

    !- End of header------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------

    flops_estielim = Real(UNAVAILABLE,KIND=8)
    flops_assemb = Real(UNAVAILABLE,KIND=8)
    flops_elim = Real(UNAVAILABLE,KIND=8)

    !---------------------------------------------------------------------------
    ! [2] Get the flops
    !---------------------------------------------------------------------------

#if HAVE_LIBMUMPS    
    If ( id_mumps%JOB >= 1 ) flops_estielim = id_mumps%rinfog(1)
    If ( id_mumps%JOB >= 2 ) flops_assemb = id_mumps%rinfog(2)
    If ( id_mumps%JOB >= 2 ) flops_elim = id_mumps%rinfog(3)
#endif

  End Subroutine XMPH_mumps_get_perfstats

  ! [+] routine : XMPH_mumps_set_schurlist -------------------------------------
  ! 
  !> Activate the schur complement computation
  !! @see XMPH_sds_set_perfstats()
  Subroutine XMPH_mumps_set_schurlist &
       (id_mumps, nschurlist , schurlist, info )

    !* Module(s) *!

    Use mph_error_mod
    Implicit None
    Include "mpif.h"

    !* Arguments *!

    Type(XMPH_mumps_t)        , Intent(inout) :: id_mumps
    MPH_INT           , Intent(in   ) :: nschurlist 
    MPH_INT           , Intent(in   ) :: schurlist (nschurlist)
    Integer              , Intent(  out) :: info

#if HAVE_LIBMUMPS    

    !* Local variables *!
    MPH_INT :: i
    Integer :: commsize

    !- End of header------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Check the parameters
    !---------------------------------------------------------------------------

    ! ncshurlist validity
    ASSRT( nschurlist >= 0 ) 

    ! current wrapper only support centralized schur
    Call MPI_Comm_size( id_mumps%COMM, commsize, info )
    ASSRT( info == MPI_SUCCESS )
    ASSRT( commsize == 1 ) 

    !---------------------------------------------------------------------------
    ! [2] Set the control paramaters
    !---------------------------------------------------------------------------

    id_mumps%ICNTL(19) = 2 ! Compute distributed schur 
    id_mumps%NPROW     = 1 
    id_mumps%NPCOL     = 1 
    id_mumps%MBLOCK    = nschurlist
    id_mumps%NBLOCK    = nschurlist

    !---------------------------------------------------------------------------
    ! [3] Set the data
    !---------------------------------------------------------------------------

    id_mumps%SIZE_SCHUR = nschurlist
    Allocate(id_mumps%LISTVAR_SCHUR(nschurlist), STAT=info)
    CHCKASSRT( info == 0, info )
    If ( info < 0) Return
    Do i=1, nschurlist
       id_mumps%LISTVAR_SCHUR(i) = schurlist(i)
    End Do

#else
    info = ERR_LIBMUMPS_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBMUMPS_UNAVAILABLE)
    Return
#endif 

  End Subroutine XMPH_mumps_set_schurlist


  ! [+] routine : XMPH_mumps_get_schur -----------------------------------------
  ! 
  !> Get the schur complement
  !! @see XMPH_sds_get_schur()
  Subroutine XMPH_mumps_get_schur (id_mumps, dmschur , info )
    

    !* Module(s) *!

    Use mph_error_mod
    Use XMPH_dense_matrix_mod
    Implicit None


    !* Arguments *!

    Type(XMPH_mumps_t)        , Intent(inout) :: id_mumps
    Type(XMPH_dense_matrix_t) , Intent(  out) :: dmschur
    Integer              , Intent(  out) :: info

#if HAVE_LIBMUMPS        

    !* Local variables *!
    Logical :: OptionSet
    Integer :: sym
    Integer :: stored

    !- End of header------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------

    ! Check if the schur computation option was set
    OptionSet= (id_mumps%ICNTL(19) == 2).OR.(id_mumps%ICNTL(19) == 3 )
    CHCKASSRT( OptionSet, info  ) 
    If ( info < 0) Return

    ! Check if we are after the factorisation
    CHCKASSRT( id_mumps%JOB >= 2, info)
    If ( info < 0) Return

    ! Check if the memory is valid
    CHCKASSRT( Associated(id_mumps%SCHUR) , info)
    If ( info < 0) Return

    ! Init output
    Call XMPH_dm_nullify(dmschur, info )
    CHCKASSRT( info >= 0, info )
    If ( info < 0) Return

    !---------------------------------------------------------------------------
    ! [2] Set the output
    !---------------------------------------------------------------------------

    ! Set the symmetry & storage flags
    Select Case ( id_mumps%SYM )
    Case (0) ! unsymmetric
       sym = DM_SYM_IsGeneral
       stored = DM_STORED_ALL
    Case (2) ! general symmetric
       sym = DM_SYM_IsSymmetric
       stored = DM_STORED_LOWER
    Case (1) ! spd
       sym = DM_SYM_IsSPD
       stored = DM_STORED_LOWER
    End Select

    ! Set the structure
    Call XMPH_dm_build(dmschur, &
         id_mumps%schur_mloc, id_mumps%schur_nloc, id_mumps%schur_lld, &
         sym,stored,id_mumps%SCHUR, info )
    CHCKASSRT( info >= 0, info )
    If ( info < 0) Return

#else
    info = ERR_LIBMUMPS_UNAVAILABLE
    Call mph_log(MSG_ERROR,LIBMUMPS_UNAVAILABLE)
    Return
#endif 

  End Subroutine XMPH_mumps_get_schur

  ! [+] routine : XMPH_mumps_free_schur -----------------------------------------
  ! 
  !> Free the schur complement
  !! @see XMPH_sds_free_schur()
  Subroutine XMPH_mumps_free_schur (id_mumps, info )

    !* Module(s) *!

    Use mph_log_mod
    Implicit None

    !* Arguments *!

    Type(XMPH_mumps_t)        , Intent(inout) :: id_mumps
    Integer              , Intent(  out) :: info

    !- End of header------------------------------------------------------------

    info = MPH_SUCCESS
    If (.Not. Associated(id_mumps%SCHUR))Then
       info = 1
       Return
    End If
    
    Deallocate(id_mumps%SCHUR)

  End Subroutine XMPH_mumps_free_schur


End module XMPH_sds_mumps_mod
