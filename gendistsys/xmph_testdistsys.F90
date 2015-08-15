! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"
! [+] program : XMPH_testdistsys ------------------------------------------------
!
!> Solve a distributed system (generated)
!! 
!!
!!----
!! Usage :
!!        <prog> <geninputfile> [<mphinputfile>] [<prefix>]
!!----
!!
Program XMPH_testdistsys

  !* Modules & co. *!
  Use mph_error_mod
  Use MPH_maphys_enum 
  Use XMPH_maphys_mod
  Use XMPH_toolkit_mod
  Use XMPH_gendistsys_mod
  Use XMPH_sparse_matrix_mod
  Implicit None
  Include 'mpif.h'

  !* Local Variable *!

  ! Parameters
  Integer, Parameter :: ndims = 3

  ! Scalars
  Integer        :: iinfo
  Integer        :: rank,commsize
  Integer        :: genunit, mphunit
  Integer        :: sym
  

  ! Arrays
  Integer        :: dims(ndims)
  XMPH_FLOAT, Pointer :: sol(:)

  ! Strings
  Character(len=MAPHYS_STRL), Parameter :: &
       FLNAME = "XMPH_ARITHmph_testdistsys.F90"

  Character(len=MAPHYS_STRL) :: geninputfile
  Character(len=MAPHYS_STRL) :: prefix
  Character(len=MAPHYS_STRL) :: mphinputfile
  Character(len=MAPHYS_STRL) :: matrixfile
  Character(len=MAPHYS_STRL) :: rhsfile
  Character(len=MAPHYS_STRL) :: outrhsfile
  Character(len=MAPHYS_STRL) :: outsolfile
  Character(len=MAPHYS_STRL) :: initguessfile

  ! Derived types
  Type(XMPH_maphys_t)        :: mphs  
  Type(XMPH_gendistsys_t)    :: gds
  Type(XMPH_sparse_matrix_t) :: sm
  
  !- End of header -------------------------------------------------------------

  !-----------------------------------------------------------------------------
  ! [-] Init
  !-----------------------------------------------------------------------------
  
  !
  iinfo = 0
  Call MPI_Init(iinfo)

  !
  Call Get_command_argument(1, geninputfile )
  Call Get_command_argument(2, mphinputfile )
  Call Get_command_argument(3, prefix )

  ! 
  mphs%comm = MPI_COMM_WORLD
  mphs%job = -1
  Call XMPH_maphys_driver(mphs)
  MPH_ONFAILURE_ABORT(mphs%iinfog(1))

  !-----------------------------------------------------------------------------
  ! [-] Read the generator parameters
  !-----------------------------------------------------------------------------

  genunit=10
  Open(UNIT=genunit, FILE=geninputfile, STATUS='old', ACTION='read')
  Read(genunit,*) gds%size_x 
  Read(genunit,*) gds%size_y
  Read(genunit,*) gds%size_z
  Read(genunit,*) gds%problem 
  Read(genunit,*) gds%an      
  Read(genunit,*) gds%p1      
  Read(genunit,*) gds%p2      
  Close(genunit)

  !
  dims = 0
  Call MPI_Comm_size(MPI_COMM_WORLD, commsize, iinfo)
  Call MPI_Comm_rank(MPI_COMM_WORLD, rank, iinfo)
  Call MPI_Dims_create(commsize,ndims,dims,iinfo)

 ! dims(1) = commsize/2
  !dims(2) = commsize/2
  !dims(3) = 1
    If (rank == 0) Write(6,*) "dims=", dims
  gds%rank      = rank
  gds%x_domains = dims(1)
  gds%y_domains = dims(2)
  gds%z_domains = dims(3)

  ASSRT(gds%x_domains*gds%y_domains*gds%z_domains == commsize )

  !-----------------------------------------------------------------------------
  ! [-] Read maphys parameters
  !-----------------------------------------------------------------------------

  ! Set the option to use distributed system on input
  mphs%icntl(ICNTL_INSYSTEM) = INSYSTEM_IsDistributed

  ! add other parameters from file
  If(Len_Trim(mphinputfile) /= 0)Then
     mphunit=11
     Open(UNIT=mphunit, FILE=mphinputfile, STATUS='old', ACTION='read')
     Call XMPH_read_param_freeformat &
          (mphunit,mphs%icntl,mphs%rcntl,sym,mphs%job,&
          matrixfile,rhsfile,initguessfile,&
          outrhsfile, outsolfile)
     Close(mphunit)


     If( Len_Trim(rhsfile   ) /= 0) rhsfile = Trim(prefix)//Trim(rhsfile)
     If (len_Trim(initguessfile) > 0) &
          initguessfile = Trim(prefix)//Trim(initguessfile)     
     If( Len_Trim(outrhsfile  ) /= 0) outrhsfile = Trim(prefix)//Trim(outrhsfile)
     If( Len_Trim(outsolfile  ) /= 0) outsolfile = Trim(prefix)//Trim(outsolfile)

  End If

  !-----------------------------------------------------------------------------
  ! [-] Set the distributed system
  !-----------------------------------------------------------------------------

  ! Get the local domain description
  Call XMPH_gendistsys_GetDomain(gds,mphs%lc_domain,mphs%comm,iinfo)
  MPH_ONFAILURE_ABORT(iinfo)

  ! Get the local matrix
  Call XMPH_gendistsys_GetMatrix (gds,sm,iinfo)
  MPH_ONFAILURE_ABORT(iinfo)

  mphs%sym    =  sm%sym
  mphs%n      =  sm%n
  mphs%nnz    =  sm%nnz
  mphs%rows   => sm%i
  mphs%cols   => sm%j
  mphs%values => sm%v

  ! Generate the second member
  ! Warning current solution is bad : no permutation at the end.
  Call XMPH_gen_sol_ones( sm%m, sol )
  Call XMPH_gen_rhs( sm, sol, mphs%rhs )
  ! export rhs to file
  Call XMPH_write_vect( mphs%rhs, mphs%n, iinfo, outrhsfile )

  ! by default perform everything
  If (mphs%job < 0) mphs%job = 6

  !-----------------------------------------------------------------------------
  ! [-] Solve it using maphys
  !-----------------------------------------------------------------------------


  Call XMPH_maphys_driver(mphs)
  MPH_ONFAILURE_ABORT(mphs%iinfog(1))
  ! Print estimation of the error commited
  If ( rank == 0 ) Write(6,*) "|A.x-b|/|b| = ", mphs%rinfog(4)

  ! export solution to file
  If (rank == 0) &
       Call XMPH_write_vect( mphs%sol, mphs%n, iinfo, outsolfile )

  ! ----------------------------------------------------------------------------
  ! [-] Exit routine
  ! ----------------------------------------------------------------------------
  
  Call MPI_Finalize(iinfo)

End Program XMPH_testdistsys


