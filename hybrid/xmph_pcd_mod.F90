! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"

#if XMPH_HAVE_ARITH_S
#define XMPH_ARITHabs Abs
#endif
#if MAPHYS_DEBUG
#define DBG(x) x
#else
#define DBG(x) ! x
#endif

! [+] module XMPH_pcd_mod ------------------------------------------------------
!
!> module to compute a preconditioner to the system.
Module XMPH_pcd_mod

  !* Module(s) *!
  Use mph_error_mod

  !* No implicit typing *!
  Implicit None

  !* Private constants *!
  Character(len=MAPHYS_STRL), Parameter, Private :: &
       FLNAME = "XMPH_ARITHmph_pcd_mod.F90"

  !* Access specifiers  *!
  Public :: XMPH_PCD_LocalExact
  Public :: XMPH_PCD_LocalApprox
  Private :: Sparsify
  Public :: XMPH_PCD_LocalApproxFromApprox


  !* Routines  *!
  Contains
    
    ! [+] routine : XMPH_PCD_LocalExact ----------------------------------------
    !
    !> Compute a local exact preconditioner
    !!
    !! Compute a local dense preconditioner "PCD" such as :
    !!
    !!  PCD = (assembled(local schur))^-1
    !! 
    !!-----
    !!
    !! @param [in,out] mphs
    !!     The maphys instance, here we access :
    !!     - dm_schur            [in]
    !!     - lc_domain           [in]
    !!     - dls_precond_schur   [out]
    !!     - iinfo/rinfo (stats) [in,out] 
    !! @param [   out] info
    !!     The routine status
    !!
    !!-----
    !! 
    !! @author Yohan Lee-tin-yien
    !! @author Azzam Haidar
    !! @author Luc   Giraud
    !!
    Subroutine XMPH_PCD_LocalExact(mphs, info )
      
      !* Module *!
      Use XMPH_maphys_type
      Use mph_error_mod
      Use MPH_mem_mod
      Use XMPH_dense_matrix_mod
      Use XMPH_dds_mod
      Use XMPH_schur_mod
      Use MPH_maphys_enum
      Implicit None

      !* Arguments *!
      Type(XMPH_maphys_t), Intent(inout) :: mphs
      Integer       , Intent(  out) :: info

      !* Local variables *!

      Real(kind=8) :: time_AsmsSchur
      Real(kind=8) :: time_PcdFacto 
      Integer(kind=8) :: mem
      Logical :: avoid_duplicate_Schur 

      !* Explicit Interface *!
      Interface
         Real(kind=8) Function MPI_WTime()
         End Function MPI_WTime
      End Interface

      ! End of header ----------------------------------------------------------

      !-------------------------------------------------------------------------
      ! [1] Set where to save the preconditioner
      !-------------------------------------------------------------------------
      !
      ! -- Optimization to reduce the memory peak --
      ! Since in we perform an implicit schur matrix vector operation.
      ! We may avoid an extra memory allocation here.
      !

      avoid_duplicate_schur = &
           ( mphs%ikeep(IKEEP_SCHUR_UselessAfter) == CURJOB_isPrecond )

      If ( avoid_duplicate_schur ) Then

         ! "mphs%dm_schur" won't be used anymore.
         ! So, we will store the preconditioner there.

         Call XMPH_dm_build( &
              mphs%dls_precond_schur%dm_A, &
              mphs%dm_schur%m,mphs%dm_schur%n,mphs%dm_schur%ld, &
              mphs%dm_schur%sym,mphs%dm_schur%stored, &
              mphs%dm_schur%v, info )
         CHCKASSRT( info >= 0, info )
         If (info < 0 ) Return

         ! access to "mphs%dm_schur" is now invalid 
         Call XMPH_dm_nullify(mphs%dm_schur, info )
         CHCKASSRT( info >= 0, info )
         If (info < 0 ) Return

      Else
         ! "mphs%dm_schur" will be used.
         ! So, we must store the preconditioner in a new memory.

         Call XMPH_dm_Dup (           & ! intents
              mphs%dm_schur,               & ! in
              mphs%dls_precond_schur%dm_A, & ! out
              info                         & ! 
              )
         CHCKASSRT( info >= 0, info )
         If (info < 0 ) Return

         Call MPH_mem_add2usage(mphs%mem,&
              XMPH_dm_sizeof(mphs%dls_precond_schur%dm_A))

      End If
      
      !------------------------------------------------------------------------- 
      ! [2] Assemble the copy of the schur
      !-------------------------------------------------------------------------

      time_AsmsSchur = MPI_WTime()
      Call XMPH_SCHUR_assemble_denseMatrix( & ! intents
           mphs%comm, mphs%lc_domain,   & ! in
           mphs%dls_precond_schur%dm_A, & ! inout
           info                         & ! out
           )
      CHCKASSRT( info >= 0, info )
      If (info < 0 ) Return
      time_AsmsSchur = MPI_WTime() - time_AsmsSchur

#if MAPHYS_DEBUG
      Call MPH_dbg_set_file("assembledschur.mtx")
      Call XMPH_dm_mmwrite(mphs%dls_precond_schur%dm_A, dbg_unit, info)
      ASSRT(info >=0)
#endif

      !-------------------------------------------------------------------------
      ! [3] Initialize dense direct solver
      !-------------------------------------------------------------------------
      
      time_pcdfacto = MPI_WTime()
      Call XMPH_dds_init (                   & ! intents
           mphs%dls_precond_schur%dds,  & ! inout
           info                         & ! out 
           )
      CHCKASSRT( info >= 0, info )
      If (info < 0 ) Return

      !-------------------------------------------------------------------------
      ! [4] Factorize the copy of the local schur 
      !-------------------------------------------------------------------------

      Call XMPH_dds_factorize(               & !
           mphs%dls_precond_schur%dds,  & ! dm_A <-- A^{-1} 
           mphs%dls_precond_schur%dm_A, & ! 
           info                         & !
           )
      CHCKASSRT( info >= 0, info )
      If (info < 0 ) Return
      time_pcdfacto = MPI_WTime() - time_pcdfacto

#if MAPHYS_DEBUG
      Call MPH_dbg_set_file("densepcd.mtx")
      Call XMPH_dm_mmwrite(mphs%dls_precond_schur%dm_A, dbg_unit, info)
      ASSRT(info >=0)
#endif

      !-------------------------------------------------------------------------
      ! [5] Save statistics
      !-------------------------------------------------------------------------

      mphs%rinfo(RINFO_TIMING_AsmsSchur) = time_AsmsSchur
      mphs%rinfo(RINFO_TIMING_PcdFacto ) = time_PcdFacto
      mphs%iinfo (IINFO_STRAT_PCDSDS  ) = -1
      mphs%iinfo (IINFO_PCD_ORDER     ) = mphs%dm_schur%n
      mphs%iinfo (IINFO_PCD_SIZEOF   ) = byte2Mbyte(mem)
      mphs%iinfo (IINFO_PCD_NBPIVOTS  ) = -1
      mphs%iinfo (IINFO_PCD_PERKEPT   ) = -1
      mphs%iinfo (IINFO_PCD_SDSMEMPEAK   ) = -1
      mphs%rinfo (RINFO_DSCHUR_RCOND1  ) = XMPH_dds_rinfo &
           (mphs%dls_precond_schur%dds, DDS_RINFO_RCOND1)
      mphs%rinfo (RINFO_DSCHUR_NORM1   ) = XMPH_dds_rinfo &
           (mphs%dls_precond_schur%dds, DDS_RINFO_NORM1)

    End Subroutine XMPH_PCD_LocalExact

    ! [+] routine : XMPH_PCD_LocalApprox ----------------------------------------
    !
    !> Compute a local approximate preconditioner
    !!
    !! Compute a local approximate preconditioner "PCD" such as :
    !!
    !!  PCD = ( sparsify( assemble(local schur)) )^-1
    !! 
    !!-----
    !!
    !! @param [in,out] mphs
    !!     The maphys instance, here we access :
    !!     - dm_schur          [in]
    !!     - lc_domain         [in]
    !!     - dls_precond_schur [out]
    !!
    !! @param [   out] info
    !!     The routine status
    !!
    !!-----
    !! 
    !! @author Yohan Lee-tin-yien
    !! @author Azzam Haidar
    !! @author Luc   Giraud
    !!
    Subroutine XMPH_PCD_LocalApprox(mphs, info )
      
      !* Module *!
      Use MPH_mem_mod
      Use MPH_maphys_enum
      Use XMPH_maphys_type
      Use XMPH_dense_matrix_mod
      Use XMPH_sparse_matrix_mod
      Use XMPH_sds_mod
      Use XMPH_sls_mod
      Use XMPH_schur_mod
      Implicit None
      Include 'mpif.h'
            

      !* Arguments *!
      Type(XMPH_maphys_t), Intent(inout) :: mphs
      Integer       , Intent(  out) :: info

      !* Local variables *!
      
      ! Scalars
      Integer :: sp_solver
      Real(kind=XMPH_FLOATKIND) :: sp_threshold
      Real(kind=8) :: time_AsmsSchur
      Real(kind=8) :: time_PcdFacto 
      Integer(kind=8) :: piv
      Real(kind=8) :: flops_estielim
      Real(kind=8) :: flops_assemb
      Real(kind=8) :: flops_elim
      Logical :: avoid_duplicate_Schur 

      ! Derived types
      Type(XMPH_dense_matrix_t        ) :: SBar
      Type(XMPH_sparse_matrix_t       ), Pointer :: SpSBar
      Type(XMPH_sds_t), Pointer :: SpSBar_1
      Type(MPH_mem_t)  :: sds_mem
      ! End of header ----------------------------------------------------------

      !-------------------------------------------------------------------------
      ! [1] Init
      !-------------------------------------------------------------------------

      sp_threshold = mphs%rkeep(RKEEP_PCD_SparsifyThreshold )
      sp_solver    = mphs%ikeep(IKEEP_SDS_Precond)

      !-------------------------------------------------------------------------
      ! [2.1] Compute SBar <- assemble (local schur)
      !-------------------------------------------------------------------------

      ! Copy
      avoid_duplicate_schur = &
           ( mphs%ikeep(IKEEP_SCHUR_UselessAfter) == CURJOB_isPrecond )

      If ( avoid_duplicate_schur ) Then

         ! "mphs%dm_schur" won't be used anymore.
         ! So, we will store the preconditioner there.

         Call XMPH_dm_build( &
              SBar, &
              mphs%dm_schur%m,mphs%dm_schur%n,mphs%dm_schur%ld, &
              mphs%dm_schur%sym,mphs%dm_schur%stored, &
              mphs%dm_schur%v, info )
         CHCKASSRT( info >= 0, info )
         If (info < 0 ) Return

         ! access to "mphs%dm_schur" is now invalid 
         Call XMPH_dm_nullify(mphs%dm_schur, info )
         CHCKASSRT( info >= 0, info )
         If (info < 0 ) Return

      Else

         Call XMPH_dm_Nullify(SBar, info)
         CHCKASSRT( info >= 0, info )
         If (info < 0 ) Return
         
         Call XMPH_dm_Dup (mphs%dm_schur, SBar, info )
         CHCKASSRT( info >= 0, info )
         If (info < 0 ) Return

         Call MPH_mem_add2usage(mphs%mem,XMPH_dm_sizeof(SBar))

      End If

      ! Assemble
      time_AsmsSchur = MPI_WTime()

      Call XMPH_SCHUR_assemble_denseMatrix( & ! intents
           mphs%comm, mphs%lc_domain,  & ! in
           SBar                     ,  & ! inout
           info                        & ! out
           )
      CHCKASSRT( info >= 0, info )
      If (info < 0 ) Return

      !-------------------------------------------------------------------------
      ! [2.2] Compute SpSBar <- sparsify(SBar)
      !-------------------------------------------------------------------------

      Allocate(SpSBar, STAT=info)
      CHCKASSRT( info == 0, info )
      If (info < 0 ) Return

      Call Sparsify(SBar, sp_threshold, SpSBar, info )
      CHCKASSRT( info == 0, info )
      If (info < 0 ) Return

      time_AsmsSchur = MPI_WTime() - time_AsmsSchur

      Call MPH_mem_add2usage(mphs%mem,XMPH_sm_sizeof(spSBar))

      ! Free memory associated to SBar
      If (avoid_duplicate_schur) Then

         Call XMPH_sds_free_schur(mphs%sls_domain%sds,info)
         MPH_ONFAILURE_RETURN(info)

      Else

         Call MPH_mem_sub2usage(mphs%mem,XMPH_dm_sizeof(SBar))
         Call XMPH_dm_Free(SBar, info)
         MPH_ONFAILURE_RETURN(info)

      End If

      !-------------------------------------------------------------------------
      ! [3] Compute PCD <- Factorize(SpSBar)
      !-------------------------------------------------------------------------
      !
      time_pcdfacto = MPI_WTime()
      !
      Allocate(SpSBar_1, STAT=info)
      CHCKASSRT( info == 0, info )
      If (info < 0 ) Return

      Call XMPH_sds_Select(sp_solver , SpSBar_1 ,info)
      CHCKASSRT( info >= 0, info )
      If (info < 0 ) Return

      Call XMPH_sds_Set_MPICommunicator(SpSBar_1,MPI_COMM_SELF,info)
      CHCKASSRT( info >= 0, info )
      If (info < 0 ) Return

      Call XMPH_sds_Set_matrix         (SpSBar_1,SpSBar,info)
      CHCKASSRT( info >= 0, info )
      If (info < 0 ) Return

      Call XMPH_sds_Analyze   (SpSBar_1,info)
      CHCKASSRT( info >= 0, info )
      If (info < 0 ) Return

      Call XMPH_sds_Factorize (SpSBar_1,info)
      CHCKASSRT( info >= 0, info )
      If (info < 0 ) Return
      
      time_pcdfacto = MPI_WTime() - time_pcdfacto

      !-------------------------------------------------------------------------
      ! [5] Save statistics
      !-------------------------------------------------------------------------

      mphs%rinfo(RINFO_TIMING_AsmsSchur) = time_AsmsSchur
      mphs%rinfo(RINFO_TIMING_PcdFacto ) = time_PcdFacto

      Call XMPH_sds_get_numstats&
           (SpSBar_1,piv)
      Call XMPH_sds_get_perfstats&
           (SpSBar_1,flops_estielim, flops_assemb, flops_elim)
      Call XMPH_sds_get_memstats&
           (SpSBar_1,sds_mem)
      Call MPH_mem_add2mem(mphs%mem,sds_mem)

      mphs%rinfo (RINFO_FLOPS_PcdEstiElim) = flops_estielim
      mphs%rinfo (RINFO_FLOPS_PcdAssemb  ) = flops_assemb
      mphs%rinfo (RINFO_FLOPS_PcdElim    ) = flops_elim
      mphs%iinfo (IINFO_STRAT_PCDSDS  ) = SpSBar_1%choice
      mphs%iinfo (IINFO_PCD_ORDER     ) = SpSBar%m

      mphs%iinfo (IINFO_PCD_SIZEOF   ) = &
           Byte2MByte(MPH_mem_getallusage(sds_mem))
      mphs%iinfo (IINFO_PCD_SDSMEMPEAK    ) = &
           Byte2MByte(MPH_mem_getallpeak   (sds_mem))

      mphs%iinfo (IINFO_PCD_NBPIVOTS  ) = INT(piv,4) 
      mphs%iinfo (IINFO_PCD_PERKEPT   ) = Ceiling( &
           100.d0*Real(SpSBar%nnz,8)/Real(mphs%dm_schur%n**2,8))

      mphs%rinfo (RINFO_DSCHUR_RCOND1  ) = -1
      mphs%rinfo (RINFO_DSCHUR_NORM1   ) = -1
      
      !-------------------------------------------------------------------------
      ! [4] Finish
      !-------------------------------------------------------------------------

      ! Save result
      Call XMPH_sls_Nullify(mphs%sls_precond_schur, info )
      CHCKASSRT( info >= 0, info )
      If (info < 0 ) Return

      mphs%sls_precond_schur%sm_A => SpSBar
      mphs%sls_precond_schur%sds  => SpSBar_1

    End Subroutine XMPH_PCD_LocalApprox

    ! [+] routine : sparsify ---------------------------------------------------
    !
    !> sparsify a dense matrix.
    !!
    !! Form the sparse matrix "sm" from the dense matrix "dm"
    !! by keeping only entries (i,j) according the values of dm, namely :
    !!  if | dm(i,j) | / | dm(i,i) | + |dm(j,j)| > threshold 
    !!
    !! The matrix "dm" must be square.
    !!
    !!----
    !!
    !! @param [in ] dm         the dense matrix
    !! @param [in ] threshold  the threshold used for sparsification
    !! @param [out] sm         the sparsified matrix
    !! @param [out] info       the routine status
    !!
    !!----
    !!
    !! @author Yohan Lee-tin-yien
    !! @author Azzam Haidar
    !! @author Luc   Giraud
    !!
    Subroutine Sparsify(dm, threshold, sm, info )

      !* Module(s) *!
      Use XMPH_dense_matrix_mod
      Use XMPH_sparse_matrix_mod
      Use mph_error_mod
      Implicit None

      !* Arguments *!
      Type(XMPH_dense_matrix_t)          , Intent(in   ) :: dm
      Real(kind=XMPH_FLOATKIND ) , Intent(in   ) :: threshold
      Type(XMPH_sparse_matrix_t)         , Intent(  out) :: sm
      Integer                       , Intent(  out) :: info

      !* Local Variables *!
      ! Scalars
      MPH_INT :: order
      MPH_INT :: nnz
      MPH_INT :: i,j,k, ij, ld
      Real(kind= XMPH_FLOATKIND) :: tmp

      ! Arrays
      Real(kind= XMPH_FLOATKIND), Pointer :: absDiag (:) ! | DM(i,i) |

      !-------------------------------------------------------------------------
      ! [1] Init variables
      !-------------------------------------------------------------------------

      ! Checks
      CHCKASSRT( dm%m >=  0, info )
      If (info < 0 ) Return
      
      CHCKASSRT( dm%m ==  dm%n, info )
      If (info < 0 ) Return

      ! only general case implemented.
      CHCKASSRT( dm%sym == DM_SYM_IsGeneral , info )
      If (info < 0 ) Return

      ! scalars
      order = dm%m
      ld    = dm%ld

      ! arrays
      Nullify(absDiag)
      Allocate( absDiag( order ), STAT= info )
      CHCKASSRT( info == 0 , info )
      If (info < 0 ) Return

      Do i=1, order
         ij = i + ld*(i-1)
         absDiag(i) = XMPH_ARITHabs( dm%v(ij) )
      End Do

      !-------------------------------------------------------------------------
      ! [2] Count the number of entries in "sm"
      !-------------------------------------------------------------------------

      nnz=0
      Do i=1,order
         Do j=1, order

            ij = i + ld*(j-1)
            tmp = XMPH_ARITHabs( dm%v(ij) ) / ( absDiag(i) + absDiag(j))
            tmp = tmp - threshold

            If ( tmp >= 0.d0 ) nnz = nnz + 1
            
         End Do
      End Do

      !-------------------------------------------------------------------------
      ! [3] Form sm
      !-------------------------------------------------------------------------

      ! create
      Call XMPH_sm_create (sm, nnz, info)
      CHCKASSRT( info >= 0 , info )
      If (info < 0 ) Goto 9999

      ! fill attributes
      sm%m   = order
      sm%n   = order
      sm%sym = dm%sym
      sm%fmt = SM_FMT_IJV

      ! fill data
      k=0
      Do i=1,order
         Do j=1, order

            ij = i + ld*(j-1)

            tmp = XMPH_ARITHabs( dm%v(ij) ) / ( absDiag(i) + absDiag(j))
            tmp = tmp - threshold

            If ( tmp >= 0.d0 ) Then 
               k=k+1
               sm%i(k) = i
               sm%j(k) = j

               sm%v(k) = DM%v(ij)
            End If

         End Do
      End Do

      ! verif
      CHCKASSRT( k == nnz , info )

      ! free
9999 Continue
      If(Associated(absDiag))Deallocate(absDiag)

    End Subroutine Sparsify


    ! [+] routine : XMPH_PCD_LocalApproxFromApprox  ----------------------------
    !
    !> Compute a local approximate preconditioner from an approximated schur
    !!
    !! Compute a local approximate preconditioner "PCD" such as :
    !!
    !!  PCD = ( assemble( approx local schur)) )^-1
    !! 
    !!-----
    !!
    !! @param [in,out] mphs
    !!     The maphys instance, here we access :
    !!     - dm_schur          [in]
    !!     - lc_domain         [in]
    !!     - dls_precond_schur [out]
    !!
    !! @param [   out] info
    !!     The routine status
    !!
    !!-----
    !! 
    !! @author Yohan Lee-tin-yien
    !! @author Azzam Haidar
    !! @author Luc   Giraud
    !!
    Subroutine XMPH_PCD_LocalApproxFromApprox(mphs, info )
     
      !* Module *!
      DBG(Use mph_dbg_mod)
      DBG(Use XMPH_sparse_matrix_mod)
      Use MPH_mem_mod
      Use MPH_maphys_enum
      Use XMPH_maphys_type
      Use XMPH_sparse_matrix_mod
      Use XMPH_sds_mod
      Use XMPH_sls_mod
      Use XMPH_schur_mod
      Use mph_log_mod
      Implicit None
      Include 'mpif.h'

      !* Arguments *!
      Type(XMPH_maphys_t), Intent(inout) :: mphs
      Integer       , Intent(  out) :: info

      !* Local variables *!
      
      ! Scalars
      Integer :: sp_solver
      Real(kind=XMPH_FLOATKIND) :: sp_threshold
      Real(kind=8) :: time_AsmsSchur
      Real(kind=8) :: time_PcdFacto 
      Integer(kind=8) :: piv
      Real(kind=8) :: flops_estielim
      Real(kind=8) :: flops_assemb
      Real(kind=8) :: flops_elim

      ! Derived types
      Type(XMPH_sparse_matrix_t), Pointer :: SpSBar
      Type(XMPH_sds_t), Pointer :: SpSBar_1
      Type(MPH_mem_t)  :: sds_mem
      ! End of header ----------------------------------------------------------

      !-------------------------------------------------------------------------
      ! [1] Init
      !-------------------------------------------------------------------------
      DBG(Call MPH_dbg_init)

      Nullify(SpSBar)
      Nullify(SpSBar_1)
      sp_threshold = mphs%rkeep(RKEEP_PCD_SparsifyThreshold )
      sp_solver    = mphs%ikeep(IKEEP_SDS_Precond)

      !-------------------------------------------------------------------------
      ! [2.1] Compute SpSBar <-- assemble (approx local schur)
      !-------------------------------------------------------------------------

      Allocate(spSBar,STAT=info)
      CHCKALLOC(info)
      If (info<0) Return

      Call XMPH_sm_nullify(spSBar,info)
      MPH_ONFAILURE_RETURN(info)

      Call XMPH_sm_dup(spSBar,mphs%sm_schur,info)
      MPH_ONFAILURE_RETURN(info)

      DBG(Call MPH_dbg_set_file("smschur2.mtx"))
      DBG(Call XMPH_sm_mmwrite(spSBar,dbg_unit,info))
      DBG(Close(dbg_unit))

      ! Assemble
      time_AsmsSchur = MPI_WTime()
      Call XMPH_SCHUR_assemble_sparseMatrix &
           (spSBar, mphs%comm, mphs%lc_domain, info )
      MPH_ONFAILURE_RETURN(info)
      time_AsmsSchur = MPI_WTime() - time_AsmsSchur

      
      DBG(Call MPH_dbg_set_file("smschur2a.mtx"))
      ! DBG(Call XMPH_sm_assemble(spSBar,spSBar%m,info))
      DBG(Call XMPH_sm_mmwrite(spSBar,dbg_unit,info))

      !-------------------------------------------------------------------------
      ! [3] Compute PCD <- Factorize(SpSBar)
      !-------------------------------------------------------------------------
      !
      time_pcdfacto = MPI_WTime()
      !
      Allocate(SpSBar_1, STAT=info)
      MPH_ONFAILURE_RETURN(info)

      Call XMPH_sds_Select(sp_solver , SpSBar_1 ,info)
      MPH_ONFAILURE_RETURN(info)

      Call XMPH_sds_Set_MPICommunicator(SpSBar_1,MPI_COMM_SELF,info)
      MPH_ONFAILURE_RETURN(info)

      Call XMPH_sds_Set_matrix         (SpSBar_1,SpSBar,info)
      MPH_ONFAILURE_RETURN(info)

      Call XMPH_sds_Analyze   (SpSBar_1,info)
      MPH_ONFAILURE_RETURN(info)

      Call XMPH_sds_Factorize (SpSBar_1,info)
      MPH_ONFAILURE_RETURN(info)

      time_pcdfacto = MPI_WTime() - time_pcdfacto

      !-------------------------------------------------------------------------
      ! [5] Save statistics
      !-------------------------------------------------------------------------

      mphs%rinfo(RINFO_TIMING_AsmsSchur) = time_AsmsSchur
      mphs%rinfo(RINFO_TIMING_PcdFacto ) = time_PcdFacto

      Call XMPH_sds_get_numstats&
           (SpSBar_1,piv)
      Call XMPH_sds_get_perfstats&
           (SpSBar_1,flops_estielim, flops_assemb, flops_elim)
      Call XMPH_sds_get_memstats&
           (SpSBar_1,sds_mem)
      Call MPH_mem_add2mem(mphs%mem,sds_mem)

      mphs%rinfo (RINFO_FLOPS_PcdEstiElim) = flops_estielim
      mphs%rinfo (RINFO_FLOPS_PcdAssemb  ) = flops_assemb
      mphs%rinfo (RINFO_FLOPS_PcdElim    ) = flops_elim
      mphs%iinfo (IINFO_STRAT_PCDSDS  ) = SpSBar_1%choice
      mphs%iinfo (IINFO_PCD_ORDER     ) = SpSBar%m

      mphs%iinfo (IINFO_PCD_SIZEOF   ) = &
           Byte2MByte(MPH_mem_getallusage(sds_mem))
      mphs%iinfo (IINFO_PCD_SDSMEMPEAK    ) = &
           Byte2MByte(MPH_mem_getallpeak   (sds_mem))

      mphs%iinfo (IINFO_PCD_NBPIVOTS  ) = INT(piv,4) 
      mphs%iinfo (IINFO_PCD_PERKEPT   ) = -1
      mphs%rinfo (RINFO_DSCHUR_RCOND1  ) = -1
      mphs%rinfo (RINFO_DSCHUR_NORM1   ) = -1
      
      !-------------------------------------------------------------------------
      ! [4] Finish
      !-------------------------------------------------------------------------

      ! Save result
      Call XMPH_sls_Nullify(mphs%sls_precond_schur, info )
      CHCKASSRT( info >= 0, info )
      If (info < 0 ) Return

      mphs%sls_precond_schur%sm_A => SpSBar
      mphs%sls_precond_schur%sds  => SpSBar_1
      DBG(Call MPH_dbg_exit)
    End Subroutine XMPH_PCD_LocalApproxFromApprox


End Module XMPH_PCD_mod
