Program testpart
  Implicit None

  ! System description
  ! ------------------
  integer x_domains,y_domains,z_domains, size_x, size_y, size_z

  Integer :: problem !> @see coorset_coef()
  Integer :: an      !> @see coorset_coef()
  Real*8  :: p1, p2, ac, ad, af  !> @see coorset_coef()

  ! MPI related
  ! -----------

  Include 'mpif.h'
  Integer :: rank

  ! The matrix
  ! -----------

  Integer :: n
  Integer :: nnz ! number of non zero
  Integer, Allocatable :: ia(:),ja(:)
  Real*8 , Allocatable :: a(:)

  ! The boundary
  ! ------------
  ! @param [out] right, left, top, bottom, up, down
  !   Flags to say which boundary is contained by the domain.
  !   0 = means no , 1 means yes
  Integer right, left, top, bottom, up, down 
  Integer, Allocatable :: iwork(:)

  ! Other
  ! -----
  Integer :: status
  Integer :: i

  ! Interface related
  ! -----------------
  Integer, Parameter :: nbViMax = 26 ! maximal number of neighbors 
  Integer :: nbVi
  Integer :: sizeIntrf
  Integer, Pointer   :: indexVi(:), ptr_Index_Intrfc(:), Index_Intrfc(:)
  Integer, Pointer   :: list_Intrfc(:)
  Integer, Pointer   :: Perm(:), invPerm(:)

  Integer :: sizeIntrf_lo
  Integer, Pointer   :: list_Intrfc_lo(:)
  Real*8, Pointer :: weight(:)

  ! External functions
  ! ------------------
  Integer, External :: getsizeInterf


  !- End of header -----------------------------

  Call MPI_Init(status)
  Call MPI_Comm_rank(MPI_COMM_WORLD, rank, status)


  !    
  !     Test Local_Discret
  !     


  x_domains = 2
  y_domains = 2
  z_domains = 2

  size_x = 3
  size_y = 3
  size_z = 3
  
  n= size_x * size_y * size_z
  nnz = 7 * n 

  Allocate ( iwork(n)) 
  sizeIntrf = getsizeInterf(rank,x_domains,y_domains,z_domains, &
       size_x,size_y,size_z,iwork)
  
  Allocate ( ia(nnz), ja(nnz), a(nnz))
  Allocate ( Perm(n), invPerm(n))

  Allocate(indexVi(nbViMax))
  Allocate(ptr_Index_Intrfc(nbViMax+1))
  Allocate(list_Intrfc(sizeIntrf))
  ! Index_Intrfc similar to list_Intrfc but references the multiplicities.
  !
  ! In Index_Intrfc :
  ! edges are referenced 3 times
  ! corners are referenced 7 times
  ! There is at most : 4 edges per direction and 8 corners shared.
  !
  Allocate(Index_Intrfc( 3*4*(size_x+size_y+size_z) + 7*8 )) 

  Allocate(list_Intrfc_lo(sizeIntrf))
  Allocate(perm(n),invPerm(n))
  Allocate(weight(sizeIntrf))

  Call setupInterf(rank,x_domains,y_domains,z_domains,size_x, size_y, &
     &              size_z,nbvi,indexVi,ptr_Index_Intrfc, Index_Intrfc, &
     &              list_Intrfc, sizeIntrf_lo, list_Intrfc_lo, &
     &              perm, invPerm, weight,iwork)

  ! 
  problem = 0 ! constant problem
  an = 0      ! isotrope
  p1 = 0.d0  ! no convection
  p2 = 0.d0  ! no convection

  Call local_Discret(rank,x_domains,y_domains,z_domains, &
       &  size_x, size_y, size_z, a, ia, ja, nnz, p1, p2,iwork, &
       &  problem, an, ac, ad, af,right, left, top, bottom,up,down)
  ! 
  Write (10+rank,*) "## domain"
  Write (10+rank,*) "# rank,right,left,top,bottom,up,down"
  Write (10+rank,*) rank,right,left,top,bottom,up,down 
  Write (10+rank,*) "## Matrix"
  Write (10+rank,*) "nnz=", nnz
  Write (10+rank,*) "# i j v"
  Do i=1,nnz
     Write(10+rank,*) ia(i), ja(i), a(i)
  End Do

  Write(*,*) "OK"

  Call MPI_Finalize()

End Program testpart

