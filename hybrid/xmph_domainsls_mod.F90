! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"

!> Handle the domain's sparse linear system
!!
!! Routines to manipulate the sparse linear system defined on each domain
!!
!! @author Yohan Lee-tin-yien
!!
Module XMPH_domainsls_mod
  
  ! Modules
  Use mph_error_mod
  Implicit none

  ! List of routines

  Public :: XMPH_DOMAINSLS_GetdmSchurAndFacto
  Public :: XMPH_DOMAINSLS_GetsmSchurAndFacto
  Public :: XMPH_DOMAINSLS_Solve_interior

  ! Private constants
  Character(len=MAPHYS_STRL), Private, Parameter :: FLNAME = &
       "XMPH_ARITHmph_domainsls_mod.F90"

Contains
  !
  ! [+] routine : XMPH_DOMAINSLS_GetdmSchurAndFacto ---------------------------------
  !> sparse direct solver factorizes the system & computes the schur complement
  !! 
  !! Compute the factors of a sparse linear system + extract the schur complement
  !!
  !!
  !! @param[in    ] comm         Controls the sparse direct solver parallelism.
  !!                             It is its MPI Communicator
  !!
  !! @param[in    ] thread_icntl Controls the sparse direct solver parallelism.
  !!                             It is its multithreading parameters.
  !!
  !! @param[in    ] which_sds    Controls which sparse direct solver to use.
  !!
  !! @param[in    ] nschur       Specifies the number of rows/columns in
  !!                             the schur complement.
  !!                             The schur complement is on the last rows/columns
  !!                             of A.
  !! 
  !! @param[in    ] A            Is the matrix of the linear system which
  !!                             bloc Aii is to be factorized. 
  !!
  !! @param[in,out] Aii_factors  Holds the factors of Aii (is the solver instance).
  !!
  !! @param[   out] dm_schur     Contains Schur complement.
  !!
  !! @param[   out] info         Specifies the routine status.  
  !!
  !! @author Azzam Haidar
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_DOMAINSLS_GetdmSchurAndFacto(    & ! intents
       comm, thread_icntl, which_sds, & ! in
       nschur, A,                     & ! 
       Aii_factors,                   & ! inout
       dm_schur,  info                & ! out
       )

    !* Module(s) *!

    Use XMPH_dense_matrix_mod
    Use XMPH_sls_mod
    Use XMPH_sds_mod
    Use mph_log_mod

    Implicit None

    !* Arguments *!

    Integer                    , Intent(in) :: comm
    Integer                    , Intent(in) :: thread_icntl(MTHREAD_ICNTL_SIZE)
    Integer                    , Intent(in) :: which_sds
    Integer                    , Intent(in) :: nschur
    Type(XMPH_sparse_matrix_t) , Intent(in) :: A
    Type(XMPH_sds_t)           , Intent(out) :: Aii_factors
    Type(XMPH_dense_matrix_t)  , Intent(out) :: dm_schur
    Integer                    , Intent(out) :: info

    !* Local Variables *!

    ! Scalars
    MPH_INT          :: nschurlist   !< size of array schurlist
    MPH_INT          :: i            !< dummy counter

    ! Arrays
    MPH_INT  , Pointer :: schurlist(:) !< list the rows/columns in the schur

    !- End of header----------------------------------------------------------

    !-------------------------------------------------------------------------
    ! [-] Initialize local variables, memory, etc.
    !-------------------------------------------------------------------------

    Nullify(schurlist)
    info = MPH_SUCCESS
    nschurlist   = nschur

    Allocate (schurlist(nschurlist), STAT = info )
    CHCKALLOC(info)
    MPH_ONFAILURE_GOTO9999(info)

    ! Here the schur are the last rows/columns of matrix.
    Do i= 1, nschurlist
       schurlist(i)  = A%m -nschur +i 
    End Do

    !-------------------------------------------------------------------------
    ! [-] Set the arguments of Aii_factors
    !-------------------------------------------------------------------------

    Call XMPH_SDS_select (which_sds,Aii_factors,info)
    MPH_ONFAILURE_GOTO9999(info)

    Call XMPH_SDS_set_MPICommunicator(Aii_factors,comm, info)
    MPH_ONFAILURE_GOTO9999(info)

    Call XMPH_SDS_set_matrix (Aii_factors,  A , info)
    MPH_ONFAILURE_GOTO9999(info)

    Call XMPH_SDS_set_schurlist(Aii_factors,nschurlist, schurlist, info) 
    MPH_ONFAILURE_GOTO9999(info)

    Call XMPH_SDS_set_multithread (Aii_factors, thread_icntl, info)
    MPH_ONFAILURE_GOTO9999(info)

    !-------------------------------------------------------------------------
    ! [-] Perform resolution until factorization
    !-------------------------------------------------------------------------

    Call XMPH_SDS_analyze   (Aii_factors,info)
    MPH_ONFAILURE_GOTO9999(info)

    Call XMPH_SDS_factorize (Aii_factors,info)
    MPH_ONFAILURE_GOTO9999(info)

    !-------------------------------------------------------------------------
    ! [-] Get the schur (data is not duplicated in memory)
    !-------------------------------------------------------------------------

    Call XMPH_SDS_get_schur(Aii_factors, dm_schur, info )
    MPH_ONFAILURE_GOTO9999(info)

    Call XMPH_dm_unsymmetrize(dm_schur)

    !-------------------------------------------------------------------------
    ! [-] Finish
    !-------------------------------------------------------------------------

9999 Continue

    ! Free memory

    If (Associated(schurlist)) Deallocate(schurlist)


  End Subroutine XMPH_DOMAINSLS_GetdmSchurAndFacto
  !
  ! [+] routine : XMPH_DOMAINSLS_GetsmSchurAndFacto --------------------------------------
  !
  !> Get a sparse Schur and the factors of the system 
  !! 
  !! Compute the factors with a sparse linear system 
  !! and the schur with pilut.
  !! The input matrix of pilut is permutated according to the one of "sds".
  !!
  !! @param[in    ] comm         Controls the sparse direct solver parallelism.
  !!                             It is its MPI Communicator.
  !! @param[in    ] thread_icntl Controls the sparse direct solver parallelism.
  !!                             It is its multithreading parameters.
  !! @param[in    ] which_sds    Controls which sparse direct solver to use.
  !! @param[in    ] Aii          Is the matrix to be factorized.
  !! @param[in,out] Aii_factors  Holds the factors of Aii (is the solver instance).
  !! @param[in    ] A            Is the matrix of the linear system which
  !!                             bloc Aii is to be factorized. 
  !! @param[in    ] nschur       Specifies the number of rows/columns in
  !!                             the schur complement.
  !!                             The schur complement is on the last rows/columns
  !!                             of A.
  !! @param[in    ] LUnzRowMax      Is the maximal number of entries in the factors.
  !! @param[in    ] LUtol        Is the tolerance in the factors.
  !! @param[in    ] SnzRowMax       Is the maximum number of entries in the schur.
  !! @param[in    ] Stol         Is the tolerance in the schur complement.
  !! @param[   out] sm_schur     The approximation of the Schur complement.
  !! @param[   out] info         Specifies the routine status.  
  !!                -1 means that pilut failed.          
  !!
  !! @author Azzam Haidar
  !! @author Yohan Lee-tin-yien
  !!
  !!
  Subroutine XMPH_DOMAINSLS_GetsmSchurAndFacto(    & 
       comm, thread_icntl, which_sds, Aii, Aii_factors, &
       nschur, A, LUnzRowMax,LUtol,SnzRowMax,Stol, & 
       sm_schur,  info )

    !* Module(s) *!
    Use mph_dbg_mod
    Use XMPH_sds_mod
    Use XMPH_sparse_matrix_mod
    Use XMPH_pilut_mod
    Implicit None

    !* Arguments *!

    Integer                    , Intent(in) :: comm
    Integer                    , Intent(in) :: thread_icntl(MTHREAD_ICNTL_SIZE)
    Integer                    , Intent(in) :: which_sds
    Type(XMPH_sparse_matrix_t) , Intent(in) :: Aii
    Type(XMPH_sds_t)           , Intent(out) :: Aii_factors

    Integer                    , Intent(in) :: nschur
    Type(XMPH_sparse_matrix_t) , Intent(in) :: A
    MPH_INT                 , Intent(in) :: LUnzRowMax
    Real(KIND=XMPH_FLOATKIND)  , Intent(in) :: LUtol
    MPH_INT                 , Intent(in) :: SnzRowMax
    Real(KIND=XMPH_FLOATKIND)  , Intent(in) :: Stol
    Type(XMPH_sparse_matrix_t) , Intent(out) :: sm_schur

    Integer                    , Intent(out) :: info

    !* Local Variables *!

    ! Scalars
    Integer :: INFO_IGNORE 
    MPH_INT :: i

    ! Array 
    MPH_INT, Pointer :: perm(:)

    ! Structures
    Type(XMPH_pilut_t) :: pilut

    !- End of header----------------------------------------------------------
    
    Nullify(perm)
    Call XMPH_pilut_init(pilut)

    !-------------------------------------------------------------------------
    ! [-] Analyze Aii
    !-------------------------------------------------------------------------

    Call XMPH_SDS_select (which_sds,Aii_factors,info)
    MPH_ONFAILURE_GOTO9999(info)

    Call XMPH_SDS_set_MPICommunicator(Aii_factors,comm, info)
    MPH_ONFAILURE_GOTO9999(info)

    Call XMPH_SDS_set_matrix (Aii_factors,  Aii , info)
    MPH_ONFAILURE_GOTO9999(info)

    Call XMPH_SDS_set_multithread (Aii_factors, thread_icntl, info)
    MPH_ONFAILURE_GOTO9999(info)    

    Call XMPH_SDS_analyze   (Aii_factors,info)
    MPH_ONFAILURE_GOTO9999(info)

    !-------------------------------------------------------------------------
    ! [-] Apply Aii permutation on A 
    !-------------------------------------------------------------------------

    Call XMPH_pilut_set_matrix(pilut, A, info )
    MPH_ONFAILURE_GOTO9999(info)    

    Allocate(perm(A%n), STAT= info )
    CHCKALLOC(info)
    If ( info < 0 ) Goto 9999 

    Call XMPH_SDS_get_permutation(Aii_factors,perm,Aii%n,info)
    MPH_ONFAILURE_GOTO9999(info)

    Do i=Aii%n+1, A%n
       perm(i) = i 
    End Do
    
    Call XMPH_pilut_set_permutation(pilut,perm,info)
    MPH_ONFAILURE_GOTO9999(info)

    !-------------------------------------------------------------------------
    ! [-] Get the schur complement 
    !-------------------------------------------------------------------------
    
    Call XMPH_pilut_set_ILULimits &
         (pilut,LUnzRowMax,LUtol,SnzRowMax,Stol,info)
    MPH_ONFAILURE_GOTO9999(info)    

    Call XMPH_pilut_set_schurorder(pilut,nschur,info)
    MPH_ONFAILURE_GOTO9999(info)

    Call XMPH_pilut_factorize(pilut,info)
    ! On failure set a specific return code
    If (info < 0 ) Then; info = -1; Goto 9999; EndIf

    Call XMPH_pilut_get_schur(pilut, sm_schur, info)
    MPH_ONFAILURE_GOTO9999(info)

    !-------------------------------------------------------------------------
    ! [-] Factorize
    !-------------------------------------------------------------------------

    Call XMPH_SDS_factorize (Aii_factors,info)
    MPH_ONFAILURE_GOTO9999(info)

    !-------------------------------------------------------------------------
    ! [-] Finish
    !-------------------------------------------------------------------------

9999 Continue
    If(Associated(perm)) Deallocate(perm)
    Call XMPH_pilut_exit(pilut)
    If (info < 0) Call XMPH_SDS_finalize( Aii_factors, INFO_IGNORE )

  End Subroutine XMPH_DOMAINSLS_GetsmSchurAndFacto


  ! [+] routine : XMPH_DOMAINSLS_Solve_interior ----------------------------------------
  !
  !> Solve the system on the local domain
  !!
  !! Solve the system on the local domain :
  !!
  !!      [Aii Aib] . [SOLi] = [RHSi]  \n
  !!      [Abi Abb]   [SOLb]   [RHSb]  \n
  !!
  !! knowing the factorisation of the local domain "Aii^{-1}",
  !! and the solution on the interface "SOLb" by using the formula :
  !! 
  !!    \f$  SOLi = Aii^{-1}.(RHSi - Aib.SOLb) \f$
  !!
  !!
  !! @param[in,out] Aii_1
  !!
  !!       The sparse direct solver instance containing the factorisation 
  !!       of  Aii.
  !!
  !! @param[in    ] Aib
  !!
  !!       The bloc Aib  of the local matrix A.
  !!
  !! @param[in    ] domain_rhs
  !!
  !!       The right-hand-side of the domain linear system.
  !!
  !! @param[in    ] bound_sol
  !!
  !!       The solution of the boundary (local schur system)
  !!
  !! @param[   out] lc_sol
  !!
  !!       The solution of the domain linear system.
  !!
  !! @param[in,out] info
  !!        
  !!       The routine status
  !!
  !!
  !! @author Azzam Haidar
  !! @author Luc   Giraud
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_DOMAINSLS_Solve_interior(     &
       Aii_1, Aib, domain_rhs, bound_sol, &
       domain_sol, info) 

    !* Module(s) used *!

    Use XMPH_sds_mod
    Use XMPH_sparse_matrix_mod
    Use XMPH_dense_matrix_mod
    Use mph_log_mod
    Implicit None
    Include 'mpif.h'

    !* Arguments *!

    Type(XMPH_sds_t), Intent(Inout) :: Aii_1
    Type(XMPH_sparse_matrix_t)       , Intent(In   ) :: Aib
    Type(XMPH_dense_matrix_t)        , Intent(In   ) :: domain_rhs
    Type(XMPH_dense_matrix_t)        , Intent(In   ) :: bound_sol
    Type(XMPH_dense_matrix_t)        , Intent(  Out) :: domain_sol
    Integer                     , Intent(  Out) :: info

    !* Local variables *!

    ! Scalars
    MPH_INT :: domain_ndof
    MPH_INT :: lim
    MPH_INT :: bound_ndof
    MPH_INT :: i  


    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [-] Initialize variables
    !---------------------------------------------------------------------------

    ! scalars
    info          = MPH_SUCCESS
    domain_ndof   = domain_rhs%m
    lim           = Aib%m
    bound_ndof    = Aib%n

    ! dense matrixes
    Call XMPH_DM_Create(domain_sol, domain_ndof, 1, domain_ndof, info )
    MPH_ONFAILURE_RETURN(info)

    !---------------------------------------------------------------------------
    ! [-] Construct the new RHS : domain_SOLi = domain_i - Aib. bound_sol
    !---------------------------------------------------------------------------
    Call XMPH_SM_VectorProduct &
         (Aib,0,-lim , bound_sol, domain_sol, info)
    MPH_ONFAILURE_RETURN(info)

    ! domain_sol%v(1:lim) = domain_rhs%v(1:lim) - domain_sol%v(1:lim)
    Do i=1,lim
       domain_sol%v(i) =  domain_rhs%v(i) - domain_sol%v(i)
    End Do

    !---------------------------------------------------------------------------
    ! [-] Solve the interior 
    !---------------------------------------------------------------------------

    Call XMPH_SDS_solve_RHS(Aii_1,        &
         domain_sol%v,             &
         domain_sol%n , domain_ndof, &
         info )
    MPH_ONFAILURE_RETURN(info)

    !---------------------------------------------------------------------------
    ! [-] Save the result ( domain_sol = (/ new_RHS Xb /)
    !---------------------------------------------------------------------------
    ! add the interface solution
    ! domain_sol%v( lim+1 : domain_ndof ) = bound_sol%v( 1 : bound_ndof )    
    Do i=1, bound_ndof
       domain_sol%v(lim+i) = bound_sol%v(i)
    End Do

  End Subroutine XMPH_DOMAINSLS_Solve_interior


End Module XMPH_domainsls_mod
