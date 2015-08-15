! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"

!> Generate a distributed system
Module XMPH_gendistsys_mod

  Use mph_error_mod
  Use MPH_part_type
  Implicit None

  ! 
  Character(len=MAPHYS_STRL), Parameter, Private :: &
       FLNAME = "XMPH_ARITHmph_gendistsys_mod.F90"
  
  
  Type XMPH_gendistsys_t
     Sequence
     ! Domain positioning 
     Integer :: rank      !< MPI rank of the domain
     Integer :: x_domains !< number of domains in direction x
     Integer :: y_domains !< number of domains in direction y
     Integer :: z_domains !< number of domains in direction z

     ! Domain description
     Integer :: size_x !< number of elements in direction x
     Integer :: size_y !< number of elements in direction y
     Integer :: size_z !< number of elements in direction z

     Integer :: right  !< domain have or not the right  boundary
     Integer :: left   !< domain have or not the left   boundary
     Integer :: top    !< domain have or not the top    boundary
     Integer :: bottom !< domain have or not the bottom boundary
     Integer :: up     !< domain have or not the up     boundary
     Integer :: down   !< domain have or not the down   boundary

     ! System description
     Integer :: problem !> @see coorset_coef()
     Integer :: an      !> @see coorset_coef()
     Real*8  :: p1, p2, ac, ad, af  !> @see coorset_coef()

     ! Workspace
     Integer, Pointer, Private :: iwork(:)

     ! Permutations to put the schur complement at the end of the matrix
     Integer, Pointer, Private :: Perm(:) 
     Integer, Pointer, Private :: invPerm(:)

  End type XMPH_gendistsys_t


  Contains


    !> Generate the domain description
    !!
    !! @param [in,out] gds the generator
    Subroutine XMPH_gendistsys_GetDomain(gds,domain,comm,info)
      
      !* Arguments *!

      Type(XMPH_gendistsys_t)    , Intent(inout) :: gds
      Type(maphys_domain_t)     , Intent(  out) :: domain
      Integer, Intent(in) :: comm
      Integer, Intent(out) :: info

      !* External functions *!
      Include 'mpif.h'
      Integer, External :: getsizeInterf
      
      !* Local variables *!

      Integer, Parameter :: nbViMax = 26 ! maximal number of neighbors 
      Integer, Parameter :: UNAVAILABLE = -1

      Integer :: n
      Integer :: nbVi
      Integer :: sizeIntrf
      Integer :: sizeIntrf_lo
      Integer :: Index_Intrfc_MaxSize

      Real*8, Pointer  :: weight(:)
      Integer, Pointer :: indexVi(:)
      Integer, Pointer :: ptr_Index_Intrfc(:)
      Integer, Pointer :: Index_Intrfc(:)
      Integer, Pointer :: list_Intrfc(:)
      Integer, Pointer :: list_Intrfc_lo(:)

      Integer :: i

      !- End of header ---------------------------------------------------------

      ! [1] Init
      ! --------------

      ! Index_Intrfc similar to list_Intrfc but references the multiplicities.
      !
      ! In Index_Intrfc :
      ! edges are referenced 3 times
      ! corners are referenced 7 times
      ! There is at most : 4 edges per direction and 8 corners shared.
      !


      n= gds%size_x * gds%size_y * gds%size_z

      Allocate ( gds%iwork(n) ) 
      sizeIntrf = getsizeInterf(&
           gds%rank,gds%x_domains,gds%y_domains,gds%z_domains, &
           gds%size_x,gds%size_y,gds%size_z,gds%iwork)

      Index_Intrfc_MaxSize = 3*4*(gds%size_x+gds%size_y+gds%size_z) + 7*8 +&
           sizeIntrf

      Allocate( gds%Perm(n), gds%invPerm(n), STAT=info )
      ASSRT(info == 0)
      Allocate(indexVi(nbViMax), ptr_Index_Intrfc(nbViMax+1), STAT=info)
      ASSRT(info == 0)
      Allocate(list_Intrfc(sizeIntrf), STAT=info)
      ASSRT(info == 0)
      Allocate(Index_Intrfc(Index_Intrfc_MaxSize), STAT=info) 
      ASSRT(info == 0)
      Allocate(list_Intrfc_lo(sizeIntrf), STAT=info)
      ASSRT(info == 0)
      Allocate(weight(sizeIntrf), STAT=info)
      ASSRT(info == 0)
       
      Call setupInterf(gds%rank,&
           & gds%x_domains,gds%y_domains,gds%z_domains,&
           & gds%size_x, gds%size_y,gds%size_z,&
           & nbvi,indexVi,ptr_Index_Intrfc, Index_Intrfc, &
           & list_Intrfc, sizeIntrf_lo, list_Intrfc_lo, &
           & gds%perm, gds%invPerm, weight,gds%iwork)

      Deallocate( list_Intrfc )

      Do i=1,ptr_Index_Intrfc(nbVi+1)-1
         Index_Intrfc(i) = gds%invPerm(Index_Intrfc(i))
      End Do

      Do i=1,sizeIntrf_lo
         list_Intrfc_lo(i) = gds%invPerm(list_Intrfc_lo(i))
      End Do

      Write(gds%rank+20,*) "perm(:)=",gds%Perm
      Write(gds%rank+20,*) "invperm(:)=",gds%invPerm

      ! Write(gds%rank+20,*) "Index_Intrfc(:) v0=",Index_Intrfc
      ! Write(gds%rank+20,*) "list_Intrfc_lo(:) v0 =",list_Intrfc_lo

      ! Fix the starting index  
      Do i=1,ptr_Index_Intrfc(nbVi+1)-1
         Index_Intrfc(i) = Index_Intrfc(i) - ( n - sizeIntrf )
      End Do

      Do i=1,sizeIntrf_lo
         list_Intrfc_lo(i) = list_Intrfc_lo(i) - ( n - sizeIntrf )
      End Do

      ! Write(gds%rank+20,*) "Index_Intrfc(:) v1 =",Index_Intrfc
      ! Write(gds%rank+20,*) "list_Intrfc_lo(:) v1=",list_Intrfc_lo

      ! [2] Set domain
      ! -----------------

      domain%myndof = n
      domain%myndofinterior = n - sizeIntrf
      domain%myndofintrf    = sizeIntrf

      domain%lenindintrf = Index_Intrfc_MaxSize
      domain%mysizeIntrf = sizeIntrf 
      domain%mynbvi      = nbVi
      domain%myindexVi => indexVi     
      domain%myptrindexVi => ptr_Index_Intrfc
      domain%myindexintrf => Index_Intrfc
      Nullify(domain%myinterface) ! list_Intrfc ! (should be removed)

      domain%myndoflogicintrf = sizeIntrf_lo
      domain%mylogicintrf => list_Intrfc_lo
      domain%weight => weight


      !==
      domain%myint_Lg    = UNAVAILABLE
      ! Call MPI_AllReduce(n,domain%gballndof,1,&
      !      MPI_INTEGER,MPI_SUM,comm,info)
      ! ASSRT(info == MPI_SUCCESS)
      Call MPI_AllReduce(sizeIntrf_lo,domain%gballintrf,1,&
           MPI_INTEGER,MPI_SUM,comm,info)
      ASSRT(info == MPI_SUCCESS)
      ! Nullify(domain%gbtoloc)     

    End Subroutine XMPH_gendistsys_GetDomain

    !> Generate the local matrix.
    !!
    !! @param[in,out] gds    the generator
    !! @param[in,out] sm     the sparse matrix.
    !! @param[   out] info  the routine status
    !!
    Subroutine XMPH_gendistsys_GetMatrix(gds, sm, info)

      !* Module(s) *!
      Use XMPH_maphys_type
      Use XMPH_sparse_matrix_mod

      !* Arguments *!


      Type(XMPH_gendistsys_t)    , Intent(inout) :: gds
      Type(XMPH_sparse_matrix_t) , Intent(out)   :: sm
      Integer, Intent(out) :: info
      
      !* Local variables *!
      Integer :: n, nnzMax, k

      !- End of header ---------------------------------------------------------

      n = gds%size_x * gds%size_y * gds%size_z
      nnzMax = n * 7               ! 7 point scheme

      Call XMPH_sm_nullify(sm, info )
      ASSRT(info == 0)

      Call XMPH_sm_create(sm, nnzMax, info )
      ASSRT(info == 0)

      sm%sym = SM_SYM_IsGeneral
      sm%m = n
      sm%n = n

      Call local_Discret(gds%rank, &
       & gds%x_domains,gds%y_domains,gds%z_domains, &
       & gds%size_x, gds%size_y,gds%size_z, &
       & sm%v, sm%i, sm%j, sm%nnz, &
       & gds%p1, gds%p2,gds%iwork,&
       & gds%problem, gds%an, gds%ac, gds%ad, gds%af, &
       & gds%right, gds%left, gds%top, gds%bottom, gds%up, gds%down)

      Do k=1,sm%nnz
         sm%i(k) = gds%invPerm(sm%i(k))
         sm%j(k) = gds%invPerm(sm%j(k))
      End Do

    End Subroutine XMPH_gendistsys_GetMatrix


  End Module XMPH_gendistsys_mod
