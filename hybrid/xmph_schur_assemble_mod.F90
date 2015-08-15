! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"

!> module to assemble the schur complement.
!!
!! 
Module XMPH_schur_assemble_mod

  Use mph_error_mod
  Implicit None

  !* Private constants *!
  Character(len=MAPHYS_STRL), Private, Parameter :: &
       FLNAME = "XMPH_schur_assemble_mod.F90"

  !* List of routines *!
  Public :: XMPH_SCHUR_Assemble_sparseMatrix

  Private :: XMPH_spSchur_getDiagBloc
  Private :: XMPH_spSchur_addDiagBloc

  Contains

   ! [+] routine : XMPH_SCHUR_assemble_sparseMatrix -----------------------------
   !
   !> Assemble the diagonal blocs of sparse local schur.
   !!
   !! Assemble the diagonal blocs of sparse local schur .
   !! (Assemble means sum the entries)
   !!
   !! @param[in,out ]  S            the sparse local schur to assemble.
   !!     - on input,  it is in IJV+CSR format without duplicates,
   !!     - on output, it is in CSR/IJV format with maybe duplicates.
   !! @param[in     ] comm          MPI communicator
   !! @param[in     ] domain        Specifies the boundaries of the domain 
   !! @param[out    ] info          the routine status
   !!
   !! @author Azzam Haidar
   !! @author Luc   Giraud
   !! @author Yohan Lee-tin-yien
   !!
   !! @version 0.2
   !!
   !! @note
   !! From version 0.2, MPI_Bsend calls was substituted with MPI_Isend ones.
   !! 
   Subroutine XMPH_SCHUR_assemble_sparseMatrix &
        (S, comm, domain, info )

    !* Module(s) *!
    Use MPH_domain_mod
    Use XMPH_sparse_matrix_mod
    Implicit None
    Include 'mpif.h'

    !* Arguments *!
    Type(XMPH_sparse_matrix_t), Intent(inout) :: S
    Integer                   , Intent(in   ) :: comm        
    Type(maphys_domain_t)     , Intent(in   ) :: domain        
    Integer                   , Intent(out  ) :: info         

    !* Local variables *!

    ! Scalars
    Integer, Parameter :: chunksize = 1024
    Integer, Parameter :: taggc = 100 ! tagg to communicate the csr component
    Integer, Parameter :: taggj = 101 ! tagg to communicate the j   component
    Integer, Parameter :: taggv = 102 ! tagg to communicate the v   component

    Integer :: i
    Integer :: neigh
    Integer :: nbNeigh
    Integer :: info_ignore
    MPH_INT :: rBlocNmax           ! maximal value of rBloc%n
    MPH_INT :: rBlocNnzAllocated
    MPH_INT :: SnnzAllocated 
    Integer :: mpistat(MPI_STATUS_SIZE)

    ! Arrays 
    MPH_INT, Pointer :: sBlocNnzAllocated(:)    ! nb entries allocated in sBlocs
    Integer, Pointer :: mpireq(:)                  ! MPI Request
    Integer, Pointer :: neighranks(:)              ! MPI Rank of neighbors

    ! Structures
    Type(XMPH_sparse_matrix_t)  :: rBloc            ! the recv diagonal bloc
    Type(XMPH_sparse_matrix_t), Pointer :: sBloc(:) ! the sent diagonal blocs

    !- End of header ---------------------------------------------------------

    !-------------------------------------------------------------------------
    ! [1] Initialize local variables
    !-------------------------------------------------------------------------

    Nullify( sBloc, sBlocNnzAllocated, mpiReq, neighranks )

    CHCKASSRT( S%fmt == SM_FMT_CSR , info )
    If( info < 0 ) Goto 9999
    !> Get the number of element allocated & fix the nnz
    !! @see XMPH_pilut_get_schur()
    SnnzAllocated = S%nnz
    S%nnz = S%csr(S%m+1)-1
    !

    Call MPH_domain_check( domain, info )
    CHCKRINFO( info )
    If ( info < 0 ) Goto 9999 

    nbNeigh = domain%mynbvi
    neighranks => domain%myIndexVi

    Allocate( sbloc(nbNeigh), sblocNnzAllocated(nbNeigh),&
         mpiReq(3*nbNeigh), STAT=info )
    CHCKALLOC( info )
    If ( info < 0 ) Goto 9999 
    
    Do neigh = 1, nbNeigh
       Call XMPH_sm_nullify(sbloc(neigh), info )
       CHCKRINFO( info )
       If ( info < 0 ) Goto 9999
    End Do

    !

    Call XMPH_sm_nullify(rbloc, info )
    CHCKRINFO( info )
    If ( info < 0 ) Goto 9999


    !-------------------------------------------------------------------------
    ! [2] Send this subdomain's contributions to each neighbor's Schur matrix 
    !-------------------------------------------------------------------------

    Do neigh =1,nbNeigh

       Call XMPH_spSchur_GetDiagBloc &
            (sbloc(neigh), sblocNnzAllocated(neigh),&
              S, domain, neigh, chunksize, info )
       CHCKRINFO( info )
       If (info < 0) Goto 9999

       ASSRT( sbloc(neigh)%nnz == sBloc(neigh)%csr(sBloc(neigh)%n+1) - 1 )

       Call MPI_ISend(sbloc(neigh)%csr,sbloc(neigh)%n+1,MPH_INTMPI, &
            neighranks(neigh),taggc, comm, mpireq(3*neigh-2), info)
       ASSRTMPI( info ) 

       Call MPI_ISend(sbloc(neigh)%j,sbloc(neigh)%nnz,MPH_INTMPI, &
            neighranks(neigh),taggj, comm, mpireq(3*neigh-1), info)
       ASSRTMPI( info ) 

       Call MPI_ISend(sbloc(neigh)%v,sbloc(neigh)%nnz,XMPH_FLOATMPI, &
            neighranks(neigh),taggv, comm, mpireq(3*neigh  ), info)
       ASSRTMPI( info ) 

    enddo

    !-------------------------------------------------------------------------
    ! [3] Receive & add neighbor's contributions to the Schur
    !-------------------------------------------------------------------------

    ! init the rbloc
    rblocNmax = 0
    rBlocNnzAllocated = 0
    Do neigh = 1, nbNeigh       
       rblocNmax= Max( MPH_domain_getBoundarySize(domain,neigh) ,rBlocNmax)
       rBlocNnzAllocated = Max( sbloc(neigh)%nnz, rBlocNnzAllocated)
    End Do
    Call XMPH_sm_create( rBloc, rBlocNnzAllocated, info )
    MPH_ONFAILURE_GOTO9999(info)

    Allocate(rBloc%csr(rblocNmax+1), STAT=info )
    CHCKALLOC(info)
    If ( info < 0 ) Goto 9999
    rBloc%fmt = SM_FMT_CSR
    rBloc%n   = 0
    rBloc%m   = 0
    rBloc%nnz = 0

    !
    Do i=1,nbNeigh
       
       Call MPI_Probe(MPI_ANY_SOURCE,taggc,comm,mpistat, info )
       ASSRTMPI( info ) 

       ! Detect which neighbor this message correspond
       Do neigh=1,nbNeigh
          If ( neighranks(neigh)  ==  mpistat(MPI_SOURCE) ) Exit
       End Do

       ! Receive the csr
       rBloc%n = MPH_domain_getBoundarySize(domain,neigh)
       Call MPI_Recv(rBloc%csr,rBloc%n+1,&
            MPH_INTMPI,neighranks(neigh),&
            taggc,comm,MPI_STATUS_IGNORE,info)
       ASSRTMPI( info )

       ! Set the nnz
       rBloc%nnz = rBloc%csr(rBloc%n+1) - 1

       ! Reallocate rBloc if there is not enough space left
       If ( rBloc%nnz > rBlocNnzAllocated ) Then

          rBloc%nnz = rBlocNnzAllocated
          rBlocNnzAllocated = (rBloc%csr(rBloc%n+1)/chunksize + 1)* chunksize
          Call XMPH_sm_realloc(rbloc,rBlocNnzAllocated,info)
          MPH_ONFAILURE_GOTO9999(info)
          rBloc%nnz = rBloc%csr(rBloc%n+1) - 1

       End If

       ! Receive the rest
       Call MPI_Recv(rBloc%j,rBloc%nnz,&
            MPH_INTMPI,neighranks(neigh),&
            taggj,comm,MPI_STATUS_IGNORE,info)
       ASSRTMPI( info )

       Call MPI_Recv(rBloc%v,rBloc%nnz,&
            XMPH_FLOATMPI,neighranks(neigh),&
            taggv,comm,MPI_STATUS_IGNORE,info)
       ASSRTMPI( info )

       ! Add the contribution

       Call XMPH_spSchur_addDiagBloc &
            (S,SnnzAllocated,rbloc,domain,neigh, chunksize, info)
       CHCKRINFO(info)
       If (info < 0) Goto 9999

    End Do

    Call MPI_Waitall(3*nbNeigh,mpireq,MPI_STATUSES_IGNORE,info)
    ASSRTMPI( info )
    info = 0
    
    !-------------------------------------------------------------------------
    ! [4] Finish
    !-------------------------------------------------------------------------
    
    If (S%fmt == SM_FMT_SPE ) Then
       Deallocate(S%csr)
       S%fmt = SM_FMT_IJV
    End If
    
9999 Continue      
    
    Nullify(neighranks)
    Call XMPH_sm_free(rBloc,info_ignore)
    Do i=1,nbNeigh
       Call XMPH_sm_free(sBloc(i),info_ignore)
    End Do
    If(Associated(sbloc)) Deallocate(sbloc, STAT=info_ignore)
    If(Associated(mpiReq)) Deallocate(mpiReq, STAT=info_ignore)
    If(Associated(sblocNnzAllocated)) &
         Deallocate(sblocNnzAllocated, STAT=info_ignore)
    
  End Subroutine XMPH_SCHUR_assemble_sparseMatrix
   !
   !----------------------------------------------------------------------------
   !
   !> Get the diagonal bloc of the schur to exchange with neighbor "neigh"
   !!
   !!
   !! @param[in,out] B     the diagonal bloc (may be reallocated)
   !! @param[in,out] BNnzAllocated   the number of allocated entries in "B"
   !! @param[in]     S     the Schur complement (sparse matrix in csr )
   !! @param[in] domain    the domain description
   !! @param[in]  neigh    the neighbor
   !! @param[in] chunksize Reallocation is done by chunk, this is its size.
   !! @param[out] info     the routine status
   !!
   Subroutine XMPH_spSchur_GetDiagBloc &
        (B, BNnzAllocated, S, domain, neigh, chunksize, info )
    
     Use MPH_domain_mod
     Use XMPH_sparse_matrix_mod
     Implicit None
     
     !* Arguments *!

     Type(XMPH_sparse_matrix_t), Intent(inout) :: B
     MPH_INT                , Intent(inout) :: BNnzAllocated
     Type(XMPH_sparse_matrix_t), Intent(in   ) :: S
     Type(MAPHYS_domain_t)        , Intent(in   ) :: domain
     Integer                   , Intent(in   ) :: neigh
     MPH_INT                , Intent(in   ) :: chunksize
     Integer                   , Intent(  out) :: info

     !* Local variables *!

     Logical     :: found
     MPH_INT  :: Sk         !< iterator on the schur
     MPH_INT  :: Bnnz       !< number of entries in the bloc
     MPH_INT  :: Brow, Bcol !< row and column indexes in the Bloc
     MPH_INT  :: Srow, Scol !< row and column indexes in the Schur

     MPH_INT, Pointer :: bound(:)
     MPH_INT          :: boundsize

     !- End of header ----------------------------------------------------------

     !--------------------------------------------------------------------------
     ! Init
     !--------------------------------------------------------------------------

     info = 0
     CHCKASSRT( S%fmt == SM_FMT_CSR , info )
     If ( info < 0 ) Goto 9999
     
     Call MPH_domain_getBoundary(domain,neigh,bound,boundsize)
     
     !--------------------------------------------------------------------------
     ! Allocate the bloc if needed
     !--------------------------------------------------------------------------

     If (B%nnz == 0)Then
        Call XMPH_sm_create(B, chunksize, info )
        MPH_ONFAILURE_GOTO9999(info)
        Bnnzallocated = B%nnz
     End If

     If ((B%fmt /= SM_FMT_CSR).OR.(B%n < boundsize))Then
        If (B%fmt == SM_FMT_CSR) Deallocate( B%csr )
        Allocate( B%csr(boundsize+1), STAT=info )
        CHCKALLOC(info)
        If ( info < 0 ) Goto 9999
        B%fmt = SM_FMT_CSR
        B%csr(1) = 1
        B%n = boundsize
        B%m = boundsize
     End If


     !--------------------------------------------------------------------------
     ! Build the bloc
     !--------------------------------------------------------------------------

     Bnnz = 0
     Do Brow =1,boundsize
        Srow = bound(Brow)

        B%csr(Brow+1) = B%csr(Brow)
        Do Bcol =1,boundsize
           Scol = bound(Bcol)

           ! Search for S(Srow,Scol)
           found = .False.
           Do Sk= S%csr(Srow),S%csr(Srow+1)-1
              If ( S%j(Sk) == Scol ) Then
                 found = .True.
                 Exit
              End If
           End Do


           ! Append S(Srow,Scol) to B
           If (found)Then

              Bnnz = Bnnz + 1

              If ( Bnnz > BnnzAllocated ) Then ! reallocate B 
                 B%nnz = Bnnz-1
                 BnnzAllocated = (Bnnz/chunksize + 1) * chunksize
                 Call XMPH_sm_realloc(B, BNnzAllocated, info )
                 MPH_ONFAILURE_GOTO9999(info)
              End If

              B%csr(Brow+1) = B%csr(Brow+1) + 1
              B%j(Bnnz) = Bcol
              B%v(Bnnz) = S%v(Sk)

           End If
        End Do
        
     End Do
     B%nnz = Bnnz

     !--------------------------------------------------------------------------
     ! Finish
     !--------------------------------------------------------------------------
           
9999 Continue     
    
     Nullify(bound)

   End Subroutine XMPH_spSchur_GetDiagBloc
   !
   !----------------------------------------------------------------------------
   !
   !> Add the diagonal bloc "B" from neighbor "neigh" to the schur "S"
   !!
   !!
   !! @param[in,out] S  the Schur complement (may be reallocated)
   !!      on input, it is in must be in csr
   !!      on output, it may be in a special file format :
   !!      it is a csr format (csr,j,v) appended
   !!      with nnz-csr(n+1) entries in coordinate formats 
   !! @param[in,out] SNnzAllocated   the number of allocated entries in "S"
   !! @param[in] B   the diagonal bloc
   !! @param[in] domain    the domain description
   !! @param[in]  neigh    the neighbor
   !! @param[in] chunksize Reallocation is done by chunk, this is its size.
   !! @param[out] info     the routine status
   !!
   Subroutine XMPH_spSchur_addDiagBloc &
        (S, SNnzAllocated, B, domain, neigh, chunksize, info )

     Use MPH_domain_mod
     Use XMPH_sparse_matrix_mod
     Implicit None
     
     !* Arguments *!

     Type(XMPH_sparse_matrix_t), Intent(inout) :: S
     MPH_INT                , Intent(inout) :: SNnzAllocated
     Type(XMPH_sparse_matrix_t), Intent(in   ) :: B
     Type(MAPHYS_domain_t)        , Intent(in   ) :: domain
     Integer                   , Intent(in   ) :: neigh
     MPH_INT                , Intent(in   ) :: chunksize
     Integer                   , Intent(  out) :: info

     !* Local variables *!

     Logical     :: found
     MPH_INT  :: Bk         !< iterator on the bloc entries
     MPH_INT  :: Sk         !< iterator on the Schur
     MPH_INT  :: Snnzin     !< number of entries in the Schur on input
     MPH_INT  :: Snnz       !< number of entries in the Schur
     MPH_INT  :: Brow, Bcol !< row and column indexes in the Bloc
     MPH_INT  :: Srow, Scol !< row and column indexes in the Schur

     MPH_INT, Pointer :: bound(:)
     MPH_INT          :: boundsize
     

     !- End of header ----------------------------------------------------------

     !--------------------------------------------------------------------------
     ! Init
     !--------------------------------------------------------------------------

     info = 0

     CHCKASSRT(B%fmt == SM_FMT_CSR, info )
     If ( info < 0 ) Goto 9999
     
     Call MPH_domain_getBoundary(domain,neigh,bound,boundsize)

     !--------------------------------------------------------------------------
     ! Add the bloc
     !--------------------------------------------------------------------------

     Snnzin = S%nnz
     Snnz = S%nnz
     Do Brow = 1,B%n

        Srow = bound(Brow)

        Do Bk = B%csr(Brow),B%csr(Brow+1)-1

           Bcol = B%j(Bk)
           Scol = bound(Bcol)

           ! Search for S(Srow,Scol)
           found = .False.
           Do Sk= S%csr(Srow),S%csr(Srow+1)-1
              If ( S%j(Sk) == Scol ) Then
                 found=.True.
                 Exit 
              End If
           End Do

           If (found)Then ! add to entry
              S%v(Sk) = S%v(Sk) + B%v(Bk) 
           Else ! create new entry at the end
              
              Snnz   = Snnz + 1
              
              If ( Snnz > SnnzAllocated ) Then ! reallocate S

                 S%nnz = Snnz-1
                 SnnzAllocated = (Snnz/chunksize + 1) * chunksize
                 Call XMPH_sm_realloc(S, SnnzAllocated, info )
                 MPH_ONFAILURE_GOTO9999(info)
                 
              End If

              S%i(Snnz) = Srow
              S%j(Snnz) = Scol
              S%v(Snnz) = B%v(Bk)

           End If

        End Do
        
     End Do

     If (Snnz /= Snnzin)Then ! S is no more in CSR
        S%nnz = Snnz
        S%fmt = SM_FMT_SPE
     End If

     !--------------------------------------------------------------------------
     ! Finish
     !--------------------------------------------------------------------------
           
9999 Continue     

     Nullify(bound)

   End Subroutine XMPH_spSchur_addDiagBloc


End Module XMPH_schur_assemble_mod
