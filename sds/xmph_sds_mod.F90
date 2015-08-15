! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"

! [+] module XMPH_sds_mod ------------------------------------------
!
!> module containing the wrappers to the sparse direct solvers
!!
!! 
Module XMPH_sds_mod

  !* Modules *!

  Use MPH_sds_enum
  Use XMPH_sparse_matrix_mod 
  Use XMPH_sds_mumps_mod
  Use XMPH_sds_pastix_mod

  !* No implicit typing *!

  Implicit None

  !* Types *!

  !> Generic structure for sparse direct solver
  Type XMPH_sds_t; Sequence
     Integer          :: choice      !< which solver to choose
     Type(XMPH_mumps_t)   :: mumps  !< mumps solver 
     Type(XMPH_pastix_t) :: pastix !< pastix solver
  End type XMPH_sds_t

  !* Access specifiers *!

  Public :: XMPH_sds_select
  Public :: XMPH_sds_set_MPIcommunicator
  Public :: XMPH_sds_set_matrix

  Public :: XMPH_sds_analyze
  Public :: XMPH_sds_factorize
  Public :: XMPH_sds_solve_rhs
  Public :: XMPH_sds_solve_sparseRHS
  Public :: XMPH_sds_finalize

  Public :: XMPH_sds_set_denseschur
  Public :: XMPH_sds_set_schurlist
  Public :: XMPH_sds_get_schur
  Public :: XMPH_sds_free_schur
  Public :: XMPH_sds_get_permutation
  Public :: XMPH_sds_set_distributed_matrix
  Public :: XMPH_sds_set_multithread
  Public :: XMPH_sds_get_numstats
  Public :: XMPH_sds_get_memstats
  Public :: XMPH_sds_get_perfstats

  !* Implementations *!

Contains
  
  ! [+] routine : XMPH_sds_select  ---------------------------------------------
  !
  !> Select which 3rd party sparse direct solver to use
  !!
  !!
  !!-----
  !!
  !! @param[in     ] selection     the 3rd party sparse direct solver to use
  !! @param[in,out ] sds           the wrapper structure
  !! @param[out    ] info          the routine status
  !!
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  !! @version 0.1 
  !!
  Subroutine XMPH_sds_select( & ! intents
       selection,       & ! in
       sds  ,           & ! inout
       info             & ! out
       )
    
    Implicit None

    !* Subroutine arguments *!

    Integer                      , Intent(in   ) :: selection
    Type(XMPH_sds_t) , Intent(inout) :: sds
    Integer                      , Intent(  out) :: info


    !- End of header------------------------------------------------------------


    Select Case(selection)
    Case( SDS_IsMUMPS ); info = XMPH_mumps_available()
    Case( SDS_IsPASTIX ); info = XMPH_pastix_available()
    Case Default ; info = -1
    End Select

    If ( info == 0 ) sds%choice = selection

  End Subroutine XMPH_sds_select

  ! [+] routine : XMPH_sds_set_MPIcommunicator ---------------------------------
  !
  !> set the MPI Communicator to be used by the direct solver
  !! 
  !! This subroutine sets
  !! the MPI communicator used by the direct solver "sds".
  !! Where the used MPI communicator is a copy of "comm".
  !!
  !!-----
  !!
  !! @param[in,out] sds  
  !!
  !!       the direct solver instance
  !!
  !! @param[in    ] comm  
  !! 
  !!       the communicator to be copied.
  !!         
  !! @param[out   ] info
  !!
  !!  the output status of the subroutine
  !!         - 0     : success 
  !!         - -1    : no direct solver available
  !!         - other : the return status of MPI_Comm_dup
  !!         .
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_sds_set_MPIcommunicator &
       (sds, comm, info )


    Implicit None
    Include 'mpif.h'
    Type(XMPH_sds_t), intent(inout) ::  sds
    Integer              , intent(in   ) :: comm
    Integer              , intent(  out) :: info

    !- End of header------------------------------------------------------------

    Select Case (sds%choice)
    Case (SDS_IsMUMPS); 
       Call XMPH_mumps_set_MPIcommunicator(sds%mumps,comm,info)
    Case (SDS_IsPASTIX);
       Call XMPH_pastix_set_MPIcommunicator(sds%pastix,comm,info)
    Case default
       info = SDS_IsUnavailable
    end select

  end subroutine XMPH_sds_set_MPIcommunicator

  ! [+] routine : XMPH_sds_Set_matrix ------------------------------------------
  !
  !> set the centralized sparse matrix to be used by the direct solver
  !!
  !! This subroutine sets
  !! the matrix of the linear system to be solved by "sds".
  !! Where the matrix is a sparse matrix centralized
  !! on the MPI process of rank 0 (according to "sds"'s MPI communicator)
  !! 
  !!-----
  !!
  !! @param[in,out] sds  
  !!
  !!       the direct solver instance
  !!
  !! @param[in    ] M  
  !!
  !!       the input matrix (only relevant on process 0)
  !!         
  !! @param[  out] info  
  !!       the output status of the subroutine
  !!         - 0     : success 
  !!         - -1    : no direct solver available,
  !!                 or another error (see below).
  !!                 see log messages to differenciate
  !!         - other : return error of the choosen 
  !!                 direct solver.
  !!                 (may contain -1, see log messages)
  !!         .
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_sds_set_matrix &
       (sds, M, info )
    Use XMPH_sparse_matrix_mod

    Implicit None
    Type(XMPH_sds_t), intent(inout) :: sds
    Type(XMPH_sparse_matrix_t), intent(in   ) :: M
    Integer              , intent(  out) :: info

    !- End of header------------------------------------------------------------

    Select Case (sds%choice)
    Case (SDS_IsMUMPS) 
       Call XMPH_mumps_set_matrix(sds%mumps,M,info)
    Case (SDS_IsPASTIX) 
       Call XMPH_pastix_set_matrix(sds%pastix,M,info)
    Case default
       info= SDS_IsUnavailable
    end select

  End Subroutine XMPH_sds_Set_Matrix


  ! [+] routine : XMPH_sds_analyze ---------------------------------------------
  !
  !> Perform the analysis step (in the resolution of the linear system)
  !!
  !!-----
  !!
  !! @param[in,out] sds    the direct solver instance
  !! @param[in,out] info  the routine status
  !!
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  subroutine XMPH_sds_analyze  &
       (sds, info )

    implicit none
    type(XMPH_sds_t), intent(inout) :: sds
    Integer              , intent(  out) :: info

    ! End of header ------------------------------------------------------------

    Select Case (sds%choice)
    Case (SDS_IsMUMPS)
       Call XMPH_mumps_analyze(sds%mumps,info)
    Case (SDS_IsPASTIX)
       Call XMPH_pastix_analyze(sds%pastix,info)
    Case default
       info= SDS_IsUnavailable
    end select

  end subroutine XMPH_sds_analyze

  ! [+] routine : XMPH_sds_factorize -------------------------------------------
  !
  !> Perform the factorization step
  !!
  !!-----
  !!
  !! @param[in,out] sds    the direct solver instance
  !! @param[in,out] info  the routine status
  !!
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  subroutine XMPH_sds_factorize (sds, info)

    implicit none
    type(XMPH_sds_t), intent(inout) :: sds
    Integer              , intent(  out) :: info

    !- End of header------------------------------------------------------------

    Select Case (sds%choice)
    Case (SDS_IsMUMPS); call XMPH_mumps_factorize(sds%mumps,info)
    Case (SDS_IsPASTIX); call XMPH_pastix_factorize(sds%pastix,info)
    Case default
       info= SDS_IsUnavailable
    end select

  end subroutine XMPH_sds_Factorize


  ! [+] routine : XMPH_sds_solve_RHS -------------------------------------------
  !
  !> perform the solve step on a dense centralized right-hand-side
  !!
  !! This subroutine launches the solve step 
  !! with the right-and-side "rhs".
  !! After this step, the solution is stored within "rhs".
  !! "rhs" is only relevant on "sds's" MPI process 0.
  !!
  !!-----
  !! 
  !! @param[in,out] sds    
  !!  
  !!        the direct solver instance
  !!
  !! @param[in,out] rhs   
  !!
  !!       - on input : the right-hand-side 
  !!                    only relevant on sds's MPI process 0.
  !!                    (in dense row major matrix format)
  !!       - on output: the solution 
  !!                    with the same comments 
  !!
  !! @param[in    ] nrhs  
  !!
  !!        number of right-hand-side to solve
  !!
  !! @param[in    ] ldrhs 
  !!
  !!        leading dimension of rhs
  !!
  !! @param[   out] info 
  !!
  !!        the output status of the subroutine
  !!        -  0    : success 
  !!        - -1    : no direct solver available,
  !!                 or another error (see below).
  !!                 see log messages to differenciate
  !!        - other : return error of the choosen 
  !!                 direct solver.
  !!                 (may contain -1, see log messages)
  !!
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_sds_solve_RHS  &
       (sds, rhs, nrhs, ldrhs, info )

    Implicit None
    Type(XMPH_sds_t), Intent(inout) :: sds
    Integer              , Intent(in   ) :: nrhs,ldrhs
    XMPH_FLOAT     , target   , Intent(inout) :: rhs(ldrhs*nrhs)
    Integer              , Intent(  out) :: info

    !- End of header -----------------------------------------------------------

    Select Case (sds%choice)
    Case (SDS_IsMUMPS) ; call XMPH_mumps_solve_RHS  &
         (sds%mumps, rhs, nrhs, ldrhs, info )
    Case (SDS_IsPASTIX); call XMPH_pastix_solve_RHS &
         (sds%pastix, rhs, nrhs, ldrhs, info )
    Case default
       info= SDS_IsUnavailable
       write(6,*) "Error : in sparse_direct_solver_solve_RHS,", &
            "no direct solver available"
    end select

  end subroutine XMPH_sds_solve_RHS

  ! [+] routine : XMPH_sds_solve_sparseRHS -------------------------------------
  !
  !> perform the solve step on a sparse centralized right-hand-side
  !!
  !! This subroutine launches the solve step 
  !! with the right-and-side "sp_rhs" (in the CSC sparse matrix format).
  !! After this step, the solution is stored
  !! in the dense row major matrix "x".
  !! "sp_rhs" and "x" are only relevant on "sds's" MPI process 0.
  !!
  !!-----
  !! 
  !! @param[in,out] sds    
  !!  
  !!        the direct solver instance
  !!
  !! @param[in    ] nx     
  !!
  !!        the number of solutions/right-hand-sides. 
  !!
  !! @param[in    ] ldx    
  !!
  !!        leading dimension of x.
  !!
  !! @param[in,out] x
  !!      
  !!          the solution. 
  !!          only relevant of "sds's" MPI process 0.
  !!          It must be allocated by the user before the subroutine call.
  !!          (in dense row major matrix format)
  !!
  !! @param[in    ] sp_rhs 
  !!
  !!          the right-hand-side in the CSC sparse matrix format.
  !!          only relevant of "sds's" MPI process 0.
  !!
  !! @param[   out] info 
  !!
  !!        the output status of the subroutine
  !!        -  0    : success 
  !!        - -1    : no direct solver available,
  !!                 or another error (see below).
  !!                 see log messages to differenciate
  !!        - other : return error of the choosen 
  !!                 direct solver.
  !!                 (may contain -1, see log messages)
  !!
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  subroutine XMPH_sds_solve_sparseRHS &
       (sds, nx, ldx, x, sp_rhs, info)
    
    Use XMPH_sparse_matrix_mod
    Implicit None
    type(XMPH_sds_t), intent(inout) :: sds
    Integer              , intent(in   ) :: nx,ldx
    XMPH_FLOAT  , target      , intent(inout) :: x(nx*ldx)

    ! warning in CSC format
    Type(XMPH_sparse_matrix_t), intent(in   ) :: sp_rhs
    Integer              , intent(  out) :: info

    !- End of header------------------------------------------------------------


    Select Case (sds%choice)
    Case (SDS_IsMUMPS) ; 
       Call XMPH_mumps_solve_sparseRHS &
            (sds%mumps, nx, ldx, x, sp_rhs, info)

    Case (SDS_IsPASTIX);
       call XMPH_pastix_solve_sparseRHS &
            (sds%pastix, nx, ldx, x, sp_rhs, info)

    Case default
       info= SDS_IsUnavailable
    end select

  end subroutine XMPH_sds_solve_sparseRHS


  ! [+] routine : XMPH_sds_finalize  -------------------------------------------
  !
  !> destroy the direct solver instance 
  !!
  !! This subroutine destroy the direct solver instance "sds".
  !! It essentially finish communications and free
  !! previously allocated memory within "sds".
  !!
  !!-----
  !!
  !! @param[in,out] sds    the direct solver instance
  !! @param[in,out] info  the routine status
  !!
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  subroutine XMPH_sds_finalize  &
       (sds, info )


    implicit none
    type(XMPH_sds_t), intent(inout) :: sds
    Integer              , intent(  out) :: info

    !- End of header------------------------------------------------------------

    Select Case (sds%choice)
    Case (SDS_IsMUMPS); call XMPH_mumps_finalize(sds%mumps,info)
    Case (SDS_IsPASTIX); call XMPH_pastix_finalize(sds%pastix,info)
    Case default
       info= SDS_IsUnavailable
    end select

  end subroutine XMPH_sds_finalize


  ! [+] routine : XMPH_sds_set_denseschur --------------------------------------
  !
  !> Ask to compute the Schur complement during the factorisation step.
  !!
  !!
  !! This a wrapper to ask a sparse direct solver 
  !! to performs the schur complement computation during the factorization step.
  !!
  !! This subroutine activate the schur computation option.
  !! It must be called before factorization ("sparse_direct_solver_factorize")
  !!
  !! The selected indexes of the matrix constituting the schur 
  !! is given by "schurlist" and the memory to put the computed
  !! schur (during factorization) is "schur" (a dense row major matrix).
  !!
  !!-----
  !! 
  !! @param[in,out] sds         the direct solver instance.
  !!                            
  !!
  !! @param[in    ] nschurlist  the size of "schurlist".\n
  !!                            It must be equal to "nschur"
  !!
  !! @param[in    ] schurlist   the list of indexes (row & column) selected 
  !!                            to constitute the schur.\n
  !!                            (size nschurlist)
  !!
  !! @param[in    ] ldschur      the leading dimension of "schur".\n
  !!                             For the moment it must be equal to "nschur".
  !!                             (see the Warning section below)
  !!                           
  !! @param[in    ] nschur       the order of "schur"
  !!
  !! @param[in,out] schur       the 1D-array that will contain
  !!                            the values of the schur matrix 
  !!                            (after factorization, by columns).\n
  !!                            It must be allocated before the subroutine call.
  !!                            (see Warning section below)
  !!
  !! @param[   out] info        the output status of the subroutine
  !!                             - 0     = success 
  !!                             - -1    = no direct solver available,
  !!                                       or another error (see below).
  !!                                       see log messages to differenciate
  !!                             - other = return error of the choosen 
  !!                                       direct solver.
  !!                                       (may contain -1, see log messages)
  !!
  !!----
  !! @warning
  !!
  !! Many direct solvers needs that "schur" is a squared matrix 
  !! of size "nschur x nschur" instead of "ldschur x nschur".
  !!
  !! The memory copy solution is not good as it is too expensive here.
  !! For the moment, the solution is to impose "ldschur == nschur".
  !!
  !!----
  !!
  !! Method
  !! - [0.0] Initialize, check arguments, etc.
  !! - [1.0] Call the selected solver
  !!
  !!----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  !! @version 0.1
  !! @date 
  !!  - Date     : Version : Comments
  !!  - 06/10/10 : 0.1     : Update the routine (Yohan Lee-tin-yien)
  !!                             - add argument checking "Use pastix/mumps_mod"
  !!                             - change coding style 
  !!                             - add debugging statements
  !!  - 2010     : 0.1     : Write the routine  (Yohan Lee-tin-yien)
  !! 
  subroutine XMPH_sds_set_denseschur(&
       sds,                     &
       schurlist, nschurlist,   &
       schur, nschur, ldschur,  &
       info                    &
       )

    !* Modules *! 

    Use XMPH_sds_mumps_mod, Only :                 &
         XMPH_mumps_set_denseschur               ! routine



    Use XMPH_sds_pastix_mod, Only :                &
         XMPH_pastix_set_denseschur              ! routine


    Implicit none

    !* Arguments *!

    type(XMPH_sds_t), intent(inout) :: sds
    Integer              , intent(in   ) :: nschurlist
    Integer              , intent(in   ) :: schurlist(nschurlist)
    Integer              , intent(in   ) :: nschur,ldschur
    XMPH_FLOAT ,pointer, intent(inout) :: schur(:)
    Integer              , intent(  out) :: info

    !- End of header------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [0] Initialize, check arguments, etc.
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Call the selected solver
    !---------------------------------------------------------------------------

    Select Case (sds%choice)

       ! mumps

    Case (SDS_IsMUMPS) ; 

       Call XMPH_MUMPS_set_denseschur               &
            (sds%mumps, schurlist, nschurlist, &
            schur, nschur, ldschur, info)
     

       ! pastix

    Case (SDS_IsPASTIX);

       Call XMPH_PASTIX_set_denseschur               &
            (sds%pastix, schurlist, nschurlist, &
            schur, nschur, ldschur, info)


       ! none available
    Case Default
       info= SDS_IsUnavailable
    End Select

  End Subroutine XMPH_sds_set_denseschur

  ! [+] routine : XMPH_sds_set_schurlist ---------------------------------------
  !
  !> Ask to compute the Schur complement during the factorisation step.
  !!
  !!
  !! This a wrapper to ask a sparse direct solver 
  !! to performs the schur complement computation during the factorization step.
  !!
  !! This subroutine activate the schur computation option.
  !! It must be called before the analysis (" XMPH_sds_analyze() ")
  !!
  !! The selected indexes of the matrix constituting the schur 
  !! is given by "schurlist". The schur complement is allocated by 
  !! the specific implementation during the factorization (" XMPH_sds_factorize() ".
  !! To obtain the schur, see XMPH_sds_get_schur().
  !!
  !!-----
  !! 
  !! @param[in,out] sds         the direct solver instance.
  !!                            
  !!
  !! @param[in    ] nschurlist  the size of "schurlist".\n
  !!                            It must be equal to "nschur"
  !!
  !! @param[in    ] schurlist   the list of indexes (row & column) selected 
  !!                            to constitute the schur.\n
  !!                            (size nschurlist)
  !!
  !! @param[   out] info        the output status of the subroutine
  !!                             - 0     = success 
  !!                             - -1    = no direct solver available,
  !!                                       or another error (see below).
  !!                                       see log messages to differenciate
  !!                             - other = return error of the choosen 
  !!                                       direct solver.
  !!                                       (may contain -1, see log messages)
  !!----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  !! @version 0.1
  !! @date 
  !!  - Date     : Version : Comments
  !!  - 14/06/11 : 0.1     : Initial implementation (Yohan Lee-tin-yien)
  !! 
  subroutine XMPH_sds_set_schurlist (&
       sds,                     &
       nschurlist, schurlist,   &
       info                    &
       )

    !* Modules *! 

    Use XMPH_sds_mumps_mod, Only :                 &
         XMPH_mumps_set_schurlist               ! routine



    Use XMPH_sds_pastix_mod, Only :                &
         XMPH_pastix_set_schurlist              ! routine


    Implicit none

    !* Arguments *!

    type(XMPH_sds_t), intent(inout) :: sds
    Integer              , intent(in   ) :: nschurlist
    Integer              , intent(in   ) :: schurlist(nschurlist)
    Integer              , intent(  out) :: info

    !- End of header------------------------------------------------------------

    Select Case (sds%choice)
    Case (SDS_IsMUMPS) ; 

       Call XMPH_MUMPS_set_schurlist               &
            (sds%mumps, nschurlist, schurlist, info)

    Case (SDS_IsPASTIX);

       Call XMPH_PASTIX_set_schurlist               &
            (sds%pastix, nschurlist, schurlist, info)

    Case Default

       info= SDS_IsUnavailable

    End Select

  End Subroutine XMPH_sds_set_schurlist

  ! [+] routine : XMPH_sds_get_schur -------------------------------------------
  !
  !> Get the computed Schur complement during the factorisation step.
  !!
  !! The schur computation option must be activated ("XMPH_sds_set_schurlist").
  !! It must be called after factorization ("XMPH_sds_factorize")
  !!
  !! The dense column major schur returned contains a pointer to the data.
  !! Thus, user must not deallocate the returned schur complement,
  !! but instead nullify it.
  !!
  !!-----
  !! 
  !! @param[in,out] sds         the direct solver instance.
  !!                            
  !!
  !! @param[in,out] dmschur     the dense schur matrix with   
  !!                dmschur%v   the 1D-array that will contain
  !!                            the values of the schur matrix 
  !!                            (after factorization, by columns).\n
  !!                            It is allocated by the solver.
  !!                            (see Warning section below)
  !!
  !! @param[   out] info        the output status of the subroutine
  !!                             - 0     = success 
  !!                             - -1    = no direct solver available,
  !!                                       or another error (see below).
  !!                                       see log messages to differenciate
  !!                             - other = return error of the choosen 
  !!                                       direct solver.
  !!                                       (may contain -1, see log messages)
  !!
  !!----
  !! @warning
  !!
  !! The solver is responsible to deallocate the schur complement matrix.
  !!
  !!----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  !! @version 0.1
  !! @date 
  !!  - Date     : Version : Comments
  !!  - 14/06/11 : 0.1     : Initial implementation (Yohan Lee-tin-yien)
  !!
  subroutine XMPH_sds_get_schur( sds, dmschur, info )

    !* Modules *! 
    Use XMPH_dense_matrix_mod, Only : &
         XMPH_dense_matrix_t      ! type
    Implicit none

    !* Arguments *!

    Type(XMPH_sds_t), Intent(inout) :: sds
    Type(XMPH_dense_matrix_t) , Intent(  out) :: dmschur
    Integer                     , Intent(  out) :: info

    !- End of header------------------------------------------------------------

    Select Case (sds%choice)
    Case (SDS_IsMUMPS ); 
       Call XMPH_MUMPS_get_schur  (sds%mumps, dmschur, info)
    Case (SDS_IsPASTIX); 
       Call XMPH_PASTIX_get_schur (sds%pastix, dmschur, info)
    Case Default
       info= SDS_IsUnavailable
    End Select

  End Subroutine XMPH_sds_get_schur

  ! [+] routine : XMPH_sds_free_schur -------------------------------------------
  !
  !> Ask the sparse direct solver to free the Schur complement
  !!
  !! - If the schur was not computed or allocated return a warning.
  !! - Access to schur obtained previously through XMPH_sds_get_schur() will be invalid.
  !!
  !!-----
  !! 
  !! @param[in,out] sds         the direct solver instance.
  !!
  !! @param[   out] info        the output status of the subroutine
  !!                             - 0     = success 
  !!                             - -1    = no direct solver available,
  !!                                       or another error (see below).
  !!                                       see log messages to differenciate
  !!                             - other = return error of the choosen 
  !!                                       direct solver.
  !!                                       (may contain -1, see log messages)
  !!
  !!----
  !! @warning
  !!
  !! The solver is responsible to deallocate the schur complement matrix.
  !!
  !!----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  !! @version 0.1
  !! @date 
  !!  - Date     : Version : Comments
  !!  - 26/10/11 : 0.1     : Initial implementation (Yohan Lee-tin-yien)
  !!
  subroutine XMPH_sds_free_schur( sds, info )

    !* Modules *! 
    Implicit none

    !* Arguments *!

    Type(XMPH_sds_t), Intent(inout) :: sds
    Integer                     , Intent(  out) :: info

    !- End of header------------------------------------------------------------

    Select Case (sds%choice)
    Case (SDS_IsMUMPS ); 
       Call XMPH_MUMPS_free_schur  (sds%mumps, info)
    Case (SDS_IsPASTIX); 
       Call XMPH_PASTIX_free_schur (sds%pastix, info)
    Case Default
       info= SDS_IsUnavailable
    End Select

  End Subroutine XMPH_sds_free_schur


  ! [+] routine : XMPH_sds_get_permutation -------------------------------------
  !
  !> Give the permutation done on the matrix by the sparse direct solver 
  !!
  !! This subroutine copy into "perm" (array of length "nperm")
  !! the permutation done by "sds" (usually done during the analysis step).
  !!
  !!-----
  !!
  !! @param[in,out] sds  
  !!
  !!       the direct solver instance
  !!
  !! @param[in    ] nperm 
  !!
  !!          the length the array "perm".
  !!         It usually correspond to the order of 
  !!         the previously given matrix.
  !!
  !! @param[in,out] perm  
  !!
  !!         the wanted permutation vector.
  !!         It must be allocated before the call to this subroutine.
  !!
  !! @param[out   ] info
  !!
  !!  the output status of the subroutine
  !!         - 0     : success 
  !!         - -1    : no direct solver available
  !!         - <0    : error
  !!         - >0    : warning
  !!         .
  !!
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  !! @version 0.1
  !! @date 
  !!  - Date     : Version : Comments
  !!  - 07/01/11 : 0.1     : Update header the routine (Yohan Lee-tin-yien)
  !!  - previous : 0.1     : Write the routine         (Yohan Lee-tin-yien)
  !!
  subroutine XMPH_sds_get_permutation &
       (sds, perm, nperm, info)

    Implicit None
    Type(XMPH_sds_t), intent(inout) ::  sds
    Integer              , intent(in   ) :: nperm
    Integer              , intent(inout) :: perm(nperm)
    Integer              , intent(  out) :: info

    !- End of header------------------------------------------------------------

    Select Case (sds%choice)
    Case (SDS_IsMUMPS) 

       Call  XMPH_mumps_get_permutation &
            (sds%mumps, perm, nperm, info)

    Case (SDS_IsPASTIX)

       Call  XMPH_pastix_get_permutation &
            (sds%pastix, perm, nperm, info)

    Case default

       info= SDS_IsUnavailable

    end select

  end subroutine XMPH_sds_get_permutation

  ! [+] routine : XMPH_sds_set_distributed_matrix ------------------------------
  !
  !> set the sparse distributed matrix of the "sds"'s linear system.
  !!
  !! This subroutine set the matrix of "sds"'s linear system,
  !! where the matrix is a sparse matrix in coordinates format,
  !! distributed throughout the "sds"'s MPI processes.
  !!
  !!----
  !!
  !! @param[in,out] sds    
  !!
  !! the direct solver instance
  !!
  !! @param[in    ] M     
  !!
  !! this MPI process part of the distributed input matrix.
  !!
  !! @param[   out] info 
  !!
  !! the output status of the subroutine
  !!         0     : success 
  !!         -1    : no direct solver available,
  !!                 or another error (see below).
  !!                 see log messages to differenciate
  !!         other : return error of the choosen 
  !!                 direct solver.
  !!                 (may contain -1, see log messages)
  !!
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  !! @version 0.1
  !! @date 
  !!  - Date     : Version : Comments
  !!  - 07/01/11 : 0.1     : Update header the routine (Yohan Lee-tin-yien)
  !!  - previous : 0.1     : Write the routine         (Yohan Lee-tin-yien)
  !!
  subroutine XMPH_sds_set_distributed_matrix &
       (sds, M, global_n, info )
    Use XMPH_sparse_matrix_mod, Only :&
         XMPH_sparse_matrix_t
    Implicit none
    Type(XMPH_sds_t), intent(inout) ::  sds
    type(XMPH_sparse_matrix_t), intent(in   ) ::  M
    Integer              , intent(in   ) :: global_n 
    Integer              , intent(  out) :: info

    !- End of header------------------------------------------------------------

    Select Case (sds%choice)
    Case (SDS_IsMUMPS) 
 
       Call XMPH_mumps_set_distributed_matrix(sds%mumps,M,global_n,info)
     
    Case (SDS_IsPASTIX)
     
       Call XMPH_pastix_set_distributed_matrix(sds%pastix,M,global_n,info)

    Case default

       info= SDS_IsUnavailable

    end select

  end subroutine XMPH_sds_set_distributed_matrix

  ! [+] routine : XMPH_sds_set_multithread -------------------------------------
  !
  !> set the multithreading 
  !!
  !! This routine activate the multi-threading facilities of 
  !! a direct solver.
  !!
  !!-----
  !!
  !! @param[in,out] sds  
  !!
  !!       the direct solver instance
  !!
  !! @param[in    ] thread_icntl(5) 
  !!
  !!    controls of multi-threadings (see below)
  !!
  !!           -   1  : binding method of threads to cores 
  !!                  - 0 : automatic binding
  !!                  - 1 : explicit  binding
  !!                  - 2 : no binding
  !!           -   2  : number of nodes
  !!           -   3  : number of cores per node
  !!           -   4  : number of threads per process
  !!           -   5  : number of processes
  !!           .
  !! @param[out   ] info
  !!
  !!    the output status of the routine
  !!
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  !! @version 0.1
  !! @date 
  !!  - Date     : Version : Comments
  !!  - 07/01/11 : 0.1     : Update header the routine (Yohan Lee-tin-yien)
  !!  - previous : 0.1     : Write the routine         (Yohan Lee-tin-yien)
  !!
  Subroutine XMPH_sds_set_multithread (sds, thread_icntl, info )

    Implicit None
    Include 'mpif.h'
    Type(XMPH_sds_t), intent(inout) :: sds
    Integer              , intent(in   ) :: thread_icntl(5)
    Integer              , intent(  out) :: info

    Integer, parameter :: WARNING_NO_MULTITHREADING = 1

    !- End of header------------------------------------------------------------
    
    Select Case (sds%choice)
    Case (SDS_IsMUMPS)

       info = WARNING_NO_MULTITHREADING

    Case (SDS_IsPASTIX)

       Call XMPH_pastix_set_multithread(sds%pastix,thread_icntl,info)

    Case default

       info = SDS_IsUnavailable

    End Select

  End Subroutine XMPH_sds_set_multithread

  ! [+] routine : XMPH_sds_get_memstats ----------------------------------------
  !
  !> get several statistics from at sparse direct solver.
  !!
  !! If a statistic is unavailable, 
  !! the correspondent component of mem should be set to "-1".
  !!
  !!-----
  !!
  !! @param[in    ] sds        the direct solver instance
  !! @param[   out] mem        the memory related of the instance.
  !!
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  !! @verbatim
  !!  - Date     : Version : Comments
  !!  - 04/05/11 : 0.1     : Write routine             (Yohan Lee-tin-yien)
  !! @endverbatim
  !!
  Subroutine XMPH_sds_get_memstats (sds, mem )

    !* Module(s) *!
    Use MPH_mem_mod, Only : &
         MPH_mem_t, & ! type
         MPH_mem_build ! routine

    !* Arguments *!

    Type(XMPH_sds_t), Intent(in ) :: sds
    Type(MPH_mem_t          ), Intent(out) :: mem

    !- End of header------------------------------------------------------------
    
    Select Case (sds%choice)
    Case (SDS_IsMUMPS) 

       Call XMPH_mumps_get_memstats(sds%mumps,mem)

    Case (SDS_IsPASTIX)

       Call XMPH_pastix_get_memstats(sds%pastix,mem)

    Case default

       Call MPH_mem_build (mem,-1_8,-1_8,-1_8,-1_8)

    End Select

  End Subroutine XMPH_sds_get_memstats



  ! [+] routine : XMPH_sds_get_numstats ----------------------------------------
  !
  !> get numerical statistics from at sparse direct solver.
  !!
  !! If a statistic is unavailable, the statistic should be set to "-1".
  !!
  !!-----
  !!
  !! @param[in    ] sds  the direct solver instance
  !! @param[   out] piv  the number of static pivots done during the factorization.
  !!                     It is only available after the call to XMPH_sds_factorize.
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  !! @verbatim
  !!  - Date     : Version : Comments
  !!  - 04/05/11 : 0.1     : Write routine             (Yohan Lee-tin-yien)
  !! @endverbatim
  !!
  Subroutine XMPH_sds_get_numstats (sds, piv )

    !* Arguments *!

    Type(XMPH_sds_t), Intent(in) :: sds
    Integer(kind=8)      , Intent(  out) :: piv

    !- End of header------------------------------------------------------------

    Select Case (sds%choice)
    Case (SDS_IsMUMPS) 

       Call XMPH_mumps_get_numstats(sds%mumps,piv)

    Case (SDS_IsPASTIX)

       Call XMPH_pastix_get_numstats(sds%pastix,piv)

    Case default

       piv = -1

    End Select

  End Subroutine XMPH_sds_get_numstats


  ! [+] routine : XMPH_sds_get_perfstats ---------------------------------------
  !
  !> get performance statistics from at sparse direct solver.
  !!
  !! If a statistic is unavailable, the statistic should be set to "-1".
  !!
  !!-----
  !!
  !! @param[in    ] sds  the direct solver instance
  !! @param[   out] flops_estielim  
  !!                The estimated floating operations for the elimination process.
  !!                It is only available after the call to XMPH_sds_analyse.
  !! @param[   out] flops_assemb  
  !!                The floating operations for the assembly process.
  !!                It is only available after the call to XMPH_sds_factorize.
  !! @param[   out] flops_elim  
  !!                The floating operations for the elimination process.
  !!                It is only available after the call to XMPH_sds_factorize.
  !!
  !!-----
  !!
  !! @author Yohan Lee-tin-yien
  !!
  !! @verbatim
  !!  - Date     : Version : Comments
  !!  - 04/05/11 : 0.1     : Write routine             (Yohan Lee-tin-yien)
  !! @endverbatim
  !!
  Subroutine XMPH_sds_get_perfstats (sds, flops_estielim, flops_assemb, flops_elim )

    !* Arguments *!

    Type(XMPH_sds_t), Intent(in) :: sds
    Real   (kind=8)      , Intent(  out) :: flops_estielim
    Real   (kind=8)      , Intent(  out) :: flops_assemb
    Real   (kind=8)      , Intent(  out) :: flops_elim

    !- End of header------------------------------------------------------------
    
    Select Case (sds%choice)
    Case (SDS_IsMUMPS) 

       Call XMPH_mumps_get_perfstats(sds%mumps,&
            flops_estielim, flops_assemb, flops_elim )

    Case (SDS_IsPASTIX)

       Call XMPH_pastix_get_perfstats(sds%pastix,&
            flops_estielim, flops_assemb, flops_elim )

    Case default
       flops_estielim = -1.d0
       flops_assemb   = -1.d0
       flops_elim     = -1.d0
    End Select

  End Subroutine XMPH_sds_get_perfstats



end module XMPH_sds_mod
