! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"

!> Auxiliary module for interface. 
!!
!!
!! @author Azzam Haidar
!! @author Luc   Giraud
!! @author Yohan Lee-tin-yien
!!
Module XMPH_schur_aux_mod
  
  Use mph_error_mod
  Implicit None

  ! List of defined routines
  Public :: XMPH_SCHUR_UpdateVector    
  Public :: XMPH_SCHUR_SchurVectorProduct
  Private :: XMPH_SCHUR_VectorProductImpl
  Public :: XMPH_SCHUR_PcdVectorProduct
  Public :: XMPH_SCHUR_VectorScalarProduct

  !* Private constants *!
  Character(len=MAPHYS_STRL), Private, Parameter :: &
       FLNAME = "XMPH_schur_aux_mod.F90"


Contains

  ! [+] routine : XMPH_SCHUR_UpdateVector -----------------------------------
  !
  !> Update the values of a distributed vector defined on the interface.
  !!
  !! The vector "x" is defined on the interface and is 
  !! distributed  into "lc_x" accross N processors.
  !!
  !! This routine simply optimize those 2 successive operations :
  !!  1-   x    <-- SUM(lc_x) 
  !!  2-   lc_x <-- distribute(x) 
  !! 
  !! @param [in,out] mphs
  !!
  !!        the structure containing :
  !!        - the interface description
  !!        - the MPI communicator and buffer used while propagating x.
  !!
  !! @param [in,out] lc_vector  
  !!
  !!        the vector to be updated.
  !!        It's size should be equal to the boundary Size.
  !!
  !! @param [   out] info
  !!
  !!        the routine status
  !!         - < 0 : Routine failed
  !!         - = 0 : Routine exited correctly
  !!         - > 0 : Routine exited correctly but with warnings
  !!         .
  !!
  !! @author Luc   Giraud
  !! @author Azzam Haidar
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_SCHUR_UpdateVector(mphs, lc_vector , info )

    !* Module(s) & co. *!
    Use XMPH_maphys_type
    Use MPH_maphys_enum, Only : RINFO_TIMING_SlvITSComm
    Use XMPH_dense_matrix_mod, Only : XMPH_dense_matrix_t
    Use mph_error_mod
    Implicit None
    include 'mpif.h'    

    !* Arguments *!
    Type(XMPH_maphys_t)       , Intent(inout) :: mphs
    Type(XMPH_dense_matrix_t) , Intent(inout) :: lc_vector
    Integer                   , Intent(  out) :: info

    !* Local Variables *!

    ! Scalars
    Integer :: tagg      !< MPI tagg
    Integer :: comm      !< MPI Communicator
    Integer :: iinfo     !< status of the called routines
    Integer :: i,k     !< dummy iterators
    Integer :: nbvi      !< field of lc_intrf. see lc_intrf%mynbvi
    Integer :: neigh     !< neighbor
    MPH_INT :: sbuffsize !< size of sbuff(:)
    MPH_INT :: rbuffsize !< size of rbuff(:)
    MPH_INT :: start 
    MPH_INT :: end
    Real(kind=8) :: timing

    ! Arrays

    Integer   , Pointer :: mpirreq(:)               !< MPI requests receivers
    Integer   , Pointer :: mpisreq(:)               !< MPI requests senders
    XMPH_FLOAT, Pointer :: sbuff (:)             !< MPI send Buffer
    XMPH_FLOAT, Pointer :: rbuff (:)             !< MPI receive Buffer
    XMPH_FLOAT, Pointer :: x    (:)              !< values of lc_vector
    Integer, Pointer    :: indexVi          (:)  !< field of boundary
    Integer, Pointer    :: ptr_Index_Intrfc (:)  !< field of boundary
    Integer, Pointer    :: Index_Intrfc     (:)  !< field of boundary

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [1] Init
    !---------------------------------------------------------------------------
    
    ! status
    info = 0
    timing = MPI_WTime()

    ! mpi
    tagg             =  44
    comm             =  mphs%comm

    ! interface description
    nbvi             =  mphs%lc_domain%mynbvi
    indexVi          => mphs%lc_domain%myIndexVi
    ptr_Index_Intrfc => mphs%lc_domain%myPtrIndexVi
    Index_Intrfc     => mphs%lc_domain%myIndexIntrf

    ! values
    x => lc_vector%v

    !---------------------------------------------------------------------------
    ! [1] Allocate the buffers
    !---------------------------------------------------------------------------
    
    Nullify( sbuff, rbuff )
    Nullify( mpirreq, mpisreq )

    Allocate( mpirreq(nbvi), mpisreq(nbvi), STAT= info )
    CHCKASSRT( info == 0, info )
    If (info < 0) Goto 9999

    sbuffsize = ptr_Index_Intrfc(nbvi+1)- ptr_Index_Intrfc(1)
    rbuffsize = ptr_Index_Intrfc(nbvi+1)- ptr_Index_Intrfc(1)
    Allocate( sbuff(sbuffsize), rbuff(rbuffsize), STAT= info )
    CHCKASSRT( info == 0, info )
    If (info < 0) Goto 9999

    Do i=1,sbuffsize 
       sbuff(i) = XMPH_FLOATZERO
       sbuff(i) = -999999999
    End Do

    Do i=1,rbuffsize 
       rbuff(i) = XMPH_FLOATZERO
    End Do

    !---------------------------------------------------------------------------
    ! [3] Start receiving 
    !---------------------------------------------------------------------------

    Do i=1,nbvi

       start = ptr_Index_Intrfc(i)
       end   = ptr_Index_Intrfc(i+1) - 1
       Call MPI_Irecv( rbuff(start),(end-start+1), &
            XMPH_FLOATMPI,indexVi(i), tagg,comm,mpirreq(i),iinfo)
       ASSRT( iinfo  == MPI_SUCCESS )

    End Do

    !---------------------------------------------------------------------------
    ! [2] Send the computed value on the interface
    !---------------------------------------------------------------------------
    
    Do i =1,nbvi

       ! Pack 
       start = ptr_Index_Intrfc(i)
       end   = ptr_Index_Intrfc(i+1) - 1
       Do k = start, end 
          sbuff(k) = x( Index_Intrfc(k) ) 
       End Do

       ! Send
       Call MPI_Isend (sbuff(start),(end-start+1),&
            XMPH_FLOATMPI,indexVi(i),tagg,comm,mpisreq(i),iinfo)
       ASSRT( iinfo  == MPI_SUCCESS )

    End Do
    
    !---------------------------------------------------------------------------
    ! [3] Unpack values from neigbhors
    !---------------------------------------------------------------------------

    ! Write(0,*) mphs%ikeep(3),"=== "
    Do i=1,nbvi
       
       Call MPI_WaitAny(nbvi, mpirreq, neigh, MPI_STATUS_IGNORE, iinfo )
       ASSRT( iinfo  == MPI_SUCCESS )

       ! Call MPI_Wait( mpirreq(i), MPI_STATUS_IGNORE, iinfo )
       ! ASSRT( iinfo  == MPI_SUCCESS )
       ! neigh = i
       ! Write(0,*) mphs%ikeep(3), "neigh = ", neigh
       start = ptr_Index_Intrfc(neigh)
       end   = ptr_Index_Intrfc(neigh+1) - 1
       Do k= start, end
          x(Index_Intrfc(k)) = x(Index_Intrfc(k)) + rbuff(k)   
       End Do

    End Do

    Call MPI_WaitAll(nbvi, mpisreq, MPI_STATUSES_IGNORE, iinfo )
    ASSRT( iinfo  == MPI_SUCCESS )

    ! --------------------------------------------------------------------------
    ! [3] Exit routine
    ! --------------------------------------------------------------------------

9999 Continue
    If (Associated(sbuff)) Deallocate(sbuff)
    If (Associated(rbuff)) Deallocate(rbuff)
    If (Associated(mpirreq)) Deallocate(mpirreq)
    If (Associated(mpisreq)) Deallocate(mpisreq)

    timing = MPI_WTime() - timing
    mphs%RINFO(RINFO_TIMING_SlvITSComm) = &
         mphs%RINFO(RINFO_TIMING_SlvITSComm) + timing

  End Subroutine XMPH_SCHUR_UpdateVector

  ! [+] routine :  XMPH_SCHUR_SchurVectorProduct ------------------------------------
  !
  !> Perform Matrix Vector Product : z <-- Schur . x
  !!
  !! Perform Matrix Vector Product : z <-- Schur . x
  !! According to the selected Schur Type.
  !! The Schur, x, and z are distributed.
  !!  
  !! @param[in,out] mphs  
  !!
  !!       the structure containing : 
  !!       - the local part of the Schur Complement (different formats)
  !!       - the data necessary to propagate the local results.
  !!
  !! @param[in,out] x      
  !!
  !!       the local part of the vector in the matrix product.
  !!
  !! @param[in,out] z      
  !!
  !!       the local part of result of the matrix product.
  !!
  !! @param[in,out] info   
  !!         
  !!       the routine status.
  !!
  !! @todo Maybe replace direct calls to blas, etc. to generic ones.
  !! @todo other strategies
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_SCHUR_SchurVectorProduct(mphs, x, z, info)

    !* Module used *!

    Use MPH_maphys_enum
    Use XMPH_maphys_type
    Use XMPH_dense_matrix_mod

    Implicit None

    !* Subroutine Arguments *!

    Type(XMPH_maphys_t)        , intent(inout) :: mphs
    Type(XMPH_dense_matrix_t)  , intent(in   ) :: x 
    Type(XMPH_dense_matrix_t)  , intent(inout) :: z 
    Integer               , intent(  out) :: info

    !* Local variables * !
    Real(kind=8) :: timing
    Real(kind=8), External :: MPI_Wtime

    !- End of header -------------------------------------------------

    timing = MPI_WTime()

    ! Perform local matrix vector product
    Select Case ( mphs%ikeep(IKEEP_ITS_MatVect) )
    Case (MAT_VECT_IsExplicit)

       Call XMPH_DM_VectorProduct(mphs%dm_schur, x, z, info )
       MPH_ONFAILURE_GOTO9999(info)

    Case (MAT_VECT_IsImplicit)

       Call XMPH_SCHUR_VectorProductImpl (&
              mphs%sm_Abi,mphs%sm_Aib,mphs%sm_Abb,&
              mphs%sls_domain%sds, x,z,info)
       MPH_ONFAILURE_GOTO9999(info)

    End Select

    ! Update the values
    Call XMPH_SCHUR_UpdateVector(mphs, z, info)
    MPH_ONFAILURE_GOTO9999(info)

    ! Finish
9999 Continue

    timing = MPI_WTime() - timing
    mphs%rinfo(RINFO_TIMING_SlvITSMatV) = timing + &
         mphs%rinfo(RINFO_TIMING_SlvITSMatV)

    ! On error, we abort since
    ! XMPH_SCHUR_UpdateVector is performs communications routine
    ASSRT( info >= 0 )

  End Subroutine XMPH_SCHUR_SchurVectorProduct

  ! [+] routine :  XMPH_SCHUR_VectorProductImpl  --------------------------------
  !
  !> Perform Schur Matrix Vector Product, but implicitely
  !!
  !! z <-- (Abi. Aii^{-1}.Aib - Abb).x
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_SCHUR_VectorProductImpl&
       (sm_Abi,sm_Aib,sm_Abb,sds_AiiF,x,z,info)

    !* Modules *!

    Use XMPH_sparse_matrix_mod
    Use XMPH_dense_matrix_mod
    Use XMPH_sds_mod
    Implicit None
    

    !* Arguments *!

    Type(XMPH_sparse_matrix_t) , Intent(in)    :: sm_Aib         
    Type(XMPH_sparse_matrix_t) , Intent(in)    :: sm_Abi         
    Type(XMPH_sparse_matrix_t) , Intent(in)    :: sm_Abb         
    Type(XMPH_sds_t)           , Intent(inout) :: sds_AiiF
    Type(XMPH_dense_matrix_t)  , Intent(in   ) :: x
    Type(XMPH_dense_matrix_t)  , Intent(inout) :: z
    Integer                    , Intent(out)   :: info 

    !* Local Variables *!

    ! scalar
    Integer    :: info_ignore
    MPH_INT :: boundSize
    MPH_INT :: ndof
    MPH_INT :: lim
    MPH_INT :: i

    ! Derived types
    Type(XMPH_dense_matrix_t) :: y  ! temporary vector
    
    !--------------------------------------------------------------------
    ! [1] Init
    !--------------------------------------------------------------------
    !> @warning about y
    !! 
    !! As the factorized system is of size ndof,
    !! y must be of size ndof.
    !!
    info = 0
    boundSize = sm_Abb%n
    ndof = sm_Aib%m + sm_Abb%n
    lim  = sm_Aib%m  
    Call XMPH_DM_Nullify(y,info_ignore)
    
    ! check that "x" is a simple vector 
    CHCKASSRT( x%n == 1, info )
    If (info < 0) Goto 9999
    
    ! y
    Call XMPH_DM_Create(y, ndof, 1, ndof, info)
    MPH_ONFAILURE_GOTO9999(info)
    
    !--------------------------------------------------------------------
    ! [2] Compute z <- Abi. Aii^{-1}.Aib.x
    !--------------------------------------------------------------------
    
    ! y <- Aib.x
    Call XMPH_SM_VectorProduct &
         ( sm_Aib, 0,-lim, x, y, info )
    MPH_ONFAILURE_GOTO9999(info)
    
    ! y <- Aii^{-1}.y
    Call XMPH_SDS_Solve_RHS               &
         ( sds_AiiF, y%v, y%n, y%ld, info )
    MPH_ONFAILURE_GOTO9999(info)
    
    ! z <- Abi.y
    Call XMPH_SM_VectorProduct&
         ( sm_Abi, -lim, 0,y, z, info )
    MPH_ONFAILURE_GOTO9999(info)
    
    !--------------------------------------------------------------------
    ! [3] Compute y <- Abb.x
    !--------------------------------------------------------------------
    
    Do i=1, y%ld*y%n
       y%v( i ) = XMPH_FLOATZERO
    End Do
    
    Call XMPH_SM_VectorProduct &
         ( sm_Abb, -lim, -lim, x, y, info )
    MPH_ONFAILURE_GOTO9999(info)
    
    !--------------------------------------------------------------------
    ! [4] Compute z <- y - z 
    !--------------------------------------------------------------------

    Do i=1,boundSize
       z%v(i) = y%v(i) - z%v(i) 
    End Do
    
    !--------------------------------------------------------------------
    ! [5] Exit 
    !--------------------------------------------------------------------

9999 Continue
    
    Call XMPH_DM_free(y,info_ignore)


  End Subroutine XMPH_SCHUR_VectorProductImpl



  ! [+] routine :  XMPH_SCHUR_PcdVectorProduct ------------------------------------
  !
  !> Apply the preconditioning  : z <-- Precond . x
  !!
  !! Perform the Matrix Vector Product, with :
  !!  - Precond : the preconditioning matrix (usually distributed)
  !!  - x       : the vector                 (usually distributed)
  !!  - z       : the result                 (usually distributed)
  !!
  !! On error this routine abort the system, to prevent MPI locking.
  !! 
  !! @param[in,out] mphs  
  !!
  !!       the structure containing : 
  !!       - the Preconditioner in different formats
  !!       - the data necessary to propagate the local results.
  !!
  !! @param[in,out] lc_x      
  !!
  !!       the (local part of the) vector in the matrix product.
  !!
  !! @param[in,out] lc_z      
  !!
  !!       the (local part of the) result of the matrix product.
  !!
  !! @param[in,out] info   
  !!         
  !!       the routine status.
  !!
  !!
  Subroutine XMPH_SCHUR_PcdVectorProduct(mphs, lc_x, lc_z, info)

    !* Module used *!
#if MAPHYS_DEBUG
    Use mph_dbg_mod
#endif
    Use XMPH_maphys_type, Only :         &
         XMPH_maphys_t                     ! type(s)
    Use XMPH_dense_matrix_mod, Only :    &
         XMPH_dense_matrix_t               ! type(s)
    Use MPH_maphys_enum
    Use XMPH_sds_mod
    Use XMPH_dds_mod
    Use mph_log_mod
    Implicit None


    !* Arguments *!

    Type(XMPH_maphys_t)        , intent(inout) :: mphs
    Type(XMPH_dense_matrix_t)  , intent(in   ) :: lc_x 
    Type(XMPH_dense_matrix_t)  , intent(inout) :: lc_z 
    Integer               , intent(  out) :: info

    !* Local Variables *!

    ! scalar
    Integer :: Pcd_Strategy
    Real(kind=8) :: timing
    Real(kind=8), External :: MPI_Wtime

    !- End of header -------------------------------------------------

    !-----------------------------------------------------------------
    ! [1] Initialize local variables
    !-----------------------------------------------------------------

    timing = MPI_WTime()
    info = 0
    Pcd_Strategy  = mphs%ikeep(IKEEP_PCD_Strategy)
    If (Pcd_Strategy == PCD_STRATEGY_isNone ) Return

    !--------------------------------------------------------------------
    ! [2] Perform Schur vector product according to the strategy
    !--------------------------------------------------------------------

    ! Copy lc_x into lc_z
    Call XMPH_ARITHcopy(lc_x%m, lc_x%v, 1,lc_z%v,1)

    ! Perform local product : lc_z <-- Pcd. lc_z

    Select Case (Pcd_Strategy)
    Case (PCD_STRATEGY_isLocalApprox, PCD_STRATEGY_isForcedByILUT) 

       Call XMPH_SDS_Solve_RHS(              &
            mphs%sls_precond_schur%sds, &
            lc_z%v,lc_z%n,lc_z%ld,    &
            info)
       CHCKRINFO(info)
       If (info < 0) Goto 9999

    Case (PCD_STRATEGY_isLocalExact) 

       Call XMPH_DDS_Solve(  &            
            mphs%dls_precond_schur%dds , & 
            mphs%dls_precond_schur%dm_A, &
            lc_z, info)
       CHCKRINFO(info)
       If (info < 0) Goto 9999

    Case Default  

       info = -1
       CHCKRINFO(info)
       If (info < 0) Goto 9999

    End Select

    !--------------------------------------------------------------------
    ! [3] Sum all lc_z entries on the interface.
    !--------------------------------------------------------------------

    Call XMPH_SCHUR_UpdateVector(mphs, lc_z, info)
    CHCKRINFO(info)
    If (info < 0) Goto 9999

    !--------------------------------------------------------------------
    ! [4] Exit routine
    !--------------------------------------------------------------------

9999 Continue

    timing = MPI_WTime() - timing
    mphs%rinfo(RINFO_TIMING_SlvITSPcdV) = timing + &
         mphs%rinfo(RINFO_TIMING_SlvITSPcdV)

    ! On error, we abort since
    ! XMPH_SCHUR_UpdateVector is performs communications routine
    ASSRT( info >= 0 )

  End Subroutine XMPH_SCHUR_PcdVectorProduct



  ! [+] routine :  XMPH_SCHUR_VectorScalarProduct ----------------------------------
  !
  !> Perform the operation  : z <-- x . y 
  !!
  !! Perform the Matrix Vector Product, with :
  !!  - x       : the first  vector          (distributed)
  !!  - y       : the second vector          (distributed)
  !!  - z       : the result                 (distributed)
  !!
  !! Where all vectors x,y,z are vectors defined on the interface.
  !!
  !! @param[in,out] mphs  
  !!
  !!       the structure containing : 
  !!       - the data necessary to propagate the local results.
  !!       - the temporary buffer on the interface "intrfbuff"
  !!       - the weight to apply on "lc_y"
  !!
  !! @param[in,out] lc_x      
  !!
  !!       the local part of "x" in the matrix product.
  !!
  !! @param[in,out] lc_y      
  !!
  !!       the local part of "y" in the matrix product.
  !!
  !! @param[in,out] lc_z      
  !!
  !!       the local part of "z" the matrix product.
  !!
  !! @param[in,out] info   
  !!         
  !!       the routine status.
  !!
  !!
  Subroutine XMPH_SCHUR_VectorScalarProduct(mphs, lc_x, lc_y, lc_z, info)

    !* Module used *!

    Use XMPH_maphys_type     , Only : &
         XMPH_maphys_t                   ! type
    Use XMPH_dense_matrix_mod, Only : &
         XMPH_dense_matrix_t,         &  ! type
         XMPH_DM_vectorproduct ! routine
    Use MPH_maphys_enum
    Implicit None
    Include "mpif.h"

    !* Subroutine Arguments *!

    Type(XMPH_maphys_t)        , intent(inout) :: mphs
    Type(XMPH_dense_matrix_t)  , intent(in   ) :: lc_x 
    Type(XMPH_dense_matrix_t)  , intent(in   ) :: lc_y
    Type(XMPH_dense_matrix_t)  , intent(inout) :: lc_z 
    Integer               , intent(  out) :: info

    !* Local Variables *!

    ! scalar
    Integer    :: iinfo
    Integer    :: comm
    Integer    :: nb_scalar_products
    MPH_INT :: i
    MPH_INT :: nrows
    Real(kind=8) :: timing

    ! Strings
    Character*1                  :: trans 

    ! Arrays
    XMPH_FLOAT, Pointer :: orthvect (:)
    XMPH_FLOAT, Pointer :: weighted_y(:)
    Real(kind=8), Pointer :: weight (:)

    !- End of header -------------------------------------------------

    !-----------------------------------------------------------------
    ! [1] Initialize local variables
    !-----------------------------------------------------------------

    info = 0
    timing = MPI_WTime()

    ! Associate input data
    comm               =  mphs%comm
    nb_scalar_products =  lc_x%n
    nrows              =  lc_y%m
    weight             => mphs%lc_domain%weight    

    ! Allocate workspace
    Nullify( orthvect, weighted_y ) 
    Allocate(&
         orthvect(nb_scalar_products),&
         weighted_y(nrows), STAT= info )
    CHCKASSRT( info == 0, info )
    If ( info < 0) Goto 9999

    !-----------------------------------------------------------------
    ! [2] Apply the weight to "lc_y" (result= weighted_y)
    !-----------------------------------------------------------------

    Do i=1,nrows
       weighted_y(i) = weight(i) * lc_y%v(i)
    End Do

    !-----------------------------------------------------------------
    ! [3] Perform the scalar product (result = orthvect)
    !-----------------------------------------------------------------
    !
    ! Here there is multiple vector/vector product to do (multiple "x")
    ! So instead of using the BLAS "xdot", we use "xmv" 
    ! with "x^T" as the matrix.
    !
    ! Temporary result is saved in orthvect
    !

    trans = 'C'
    Call XMPH_ARITHgemv  (   & !                              
         trans,lc_x%m,       & !   orthvect(i) = x(:,i)^T . y 
         nb_scalar_products, & !                              
         XMPH_FLOATONE,      & ! 
         lc_x%v(1),lc_x%ld,  & !                              
         weighted_y(1) ,1,   & ! With                         
         XMPH_FLOATZERO,     & !   i = 1 .. nb_scalar_products
         orthvect(1), 1      & !
         )

    !-----------------------------------------------------------------
    ! [4] Compute "lc_z" , sum of "orthvect" on all processor
    !-----------------------------------------------------------------

    Call MPI_AllReduce(      & ! 
         orthvect,lc_z%v,    & ! lc_z(i) = SUM( orthvect(i) ) 
         nb_scalar_products, & !            
         XMPH_FLOATMPI,      & ! With 
         MPI_SUM,comm,iinfo  & !  i = 1 .. nb_scalar_products
         )
    ASSRT( iinfo == MPI_SUCCESS )

    !-----------------------------------------------------------------
    ! [5] Exit routine
    !-----------------------------------------------------------------

9999 Continue

    ! Free buffers
    If (Associated(orthvect)) Deallocate(orthvect)
    If (Associated(weighted_y)) Deallocate(weighted_y)

    timing = MPI_WTime() - timing
    mphs%rinfo(RINFO_TIMING_SlvITSDotP) = timing + &
         mphs%rinfo(RINFO_TIMING_SlvITSDotP)

  End Subroutine XMPH_SCHUR_VectorScalarProduct


End Module XMPH_SCHUR_aux_mod




