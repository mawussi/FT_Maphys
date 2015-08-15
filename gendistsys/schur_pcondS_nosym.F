       program scalable_mumps
*
       implicit none
*
       real*8 , dimension (:), allocatable, target :: a 
       real*8 , dimension (:), allocatable :: weight, buff
       integer, dimension (:), allocatable, target ::  ia, ja
       integer, dimension (:), allocatable :: iwork

       integer, dimension (:), allocatable :: iabis, jabis
       real*8,  dimension (:), allocatable :: abis
       real*8 , dimension (:), allocatable :: wor
       integer, dimension (:), allocatable :: iwor
       real*8   Rcond,Ferr,Berr,Rcondmin,Rcondmax
       real*8   Normm
       integer  infolap
       real*8   Normmsol,Normmbis,Norm_intrf,Norm_int

*
       integer  nbdomain
       integer, dimension (:)  , allocatable ::  IPIV  
       real*8 , dimension (:  ), allocatable :: tempAcR,tempAcRbis 
       real*8 , dimension (:,:), allocatable :: Acoa, Rt
*
       integer  scaling
       integer, dimension (:)  , allocatable ::  IPIV_S  
       real*8 , dimension (:), allocatable, target ::  S 
       real*8, dimension (:), allocatable :: SBar,sqdiag_S
*
       integer nbviM
       parameter (nbviM=27)
       integer  ptr_Index_Intrfc(nbviM+1), indexVi(nbviM+1)
       integer  sizeIntrf, sizeIntrf_lo
       integer, dimension (:), allocatable :: Index_Intrfc, list_Intrfc
       integer, dimension (:), allocatable :: Perm, invPerm ,per
       integer, dimension (:), allocatable :: list_Intrfc_lo
*
       integer nnzabb, nnzS
       integer, dimension (:), allocatable :: iabb, jabb, is, js
       real*8,  dimension (:), allocatable ::  abb, sii 
*
       integer nnzabi, nnzaib
       integer, dimension (:), allocatable :: iabi, jabi
       integer, dimension (:), allocatable :: iaib, jaib
       real*8,  dimension (:), allocatable :: aib, abi
*
       integer Up, Down, Top, Bottom, Right, Left
       integer size_x, size_y, size_z, me,x_domains,y_domains,z_domains
       integer nbvi, nnz, i, j, n, nxy, nglob, k,l 
*
       real*8, dimension (:), allocatable :: x, yaux, workaux, xaux,xsol
       real*8, dimension (:), allocatable :: orthvect
       real*8, dimension (:), allocatable, target :: temp, b
       integer lwork
       integer problem, an
       real*8  p1,p2,ac,ad,af
*
       real*8  ddotw, result, res(128)
       real*8  ZERO, ONE
       parameter (ZERO = 0.0d0, ONE = 1.0d0)
       integer nps, nmv, npr
       integer irc(5), icntl(8), nbscal,restrt,ihist, info(8) 
       integer colx, coly, colz, iconf(8), revcom, ios,ierr
       real*8  cntl(5), rinfo(4), conf(6)
       real*8  somme
       character*20 fileout
*
* Timing variable
*
       real*8  avertot,mintot,maxtot
       real*8, dimension(:), allocatable :: o_tot,o_tot_bis

       real*8  tcg, tpreloc, tfacto, tsymb, tsolv,tsolvi, taux, ttotal
       real*8  tassemb, tprefact, tdot, tMatVect, tpreapp, ttot,tset

       real*8  avercg, averpreloc, averfact, aversymb, aversolv,avertime
       real*8  averassemb, averprefact, averMatVect, averpreapp,averdot

       real*8  mincg, minpreloc, minfact, minsymb, minsolv,mintime
       real*8  minassemb, minprefact, minMatVect, minpreapp,mindot

       real*8  maxcg, maxpreloc, maxfact, maxsymb, maxsolv,maxtime
       real*8  maxassemb, maxprefact, maxMatVect, maxpreapp,maxdot

       integer minmemo,maxmemo,avermemo
       integer mineffmemo,maxeffmemo,avereffmemo
       integer minschur,maxschur,averschur
       integer, dimension(:), allocatable :: o_memo,o_memo_bis 
       integer, dimension(:), allocatable :: o_effmemo,o_effmemo_bis 
       integer, dimension(:), allocatable :: o_schur,o_schur_bis

       real*8, dimension(:), allocatable :: o_tsymb,o_tsymb_bis 
       real*8, dimension(:), allocatable :: o_tfacto,o_tfacto_bis
       real*8, dimension(:), allocatable :: o_tsolv,o_tsolv_bis
       real*8, dimension(:), allocatable :: o_tpreloc,o_tpreloc_bis
       real*8, dimension(:), allocatable :: o_tassemb,o_tassemb_bis
       real*8, dimension(:), allocatable :: o_tprefact,o_tprefact_bis
       real*8, dimension(:), allocatable :: o_tcg,o_tcg_bis
       real*8, dimension(:), allocatable :: o_tMatVect,o_tMatVect_bis
       real*8, dimension(:), allocatable :: o_tpreapp,o_tpreapp_bis
       real*8, dimension(:), allocatable :: o_tdot,o_tdot_bis
 
*
* communication buffer
*
       real*8, dimension(:), allocatable :: combuff
*
* fileoutput variable
*
       integer myIntrf(2), maxminIntrf(2,2),Job_id
*
* MUMPS data structure
*
       INCLUDE 'dmumps_struc.h'
       TYPE (DMUMPS_STRUC) mumps_seq
*
       real*8   ddotwS
       external ddotwS
*
       real*8   DDOT
       external DDOT
*
       integer getsizeInterf
       external getsizeInterf
*
C       real*8   dlange, dlansy
C       external dlange, dlansy
*
       include 'mpif.h'
       integer max_neigintrf
       integer totalInterf,totalInt
       integer nproc, infompi, comm, sizebuff, totintrf,totintrfbis
       intrinsic sqrt,max
*
*
      ttot     = ZERO
      ttotal   = ZERO
      tpreapp  = ZERO
      tdot     = ZERO
      tMatVect = ZERO
      tassemb  = ZERO
      tprefact = ZERO
      tpreloc  = ZERO
      tcg      = ZERO
      tfacto   = ZERO
      tsolv    = ZERO
      tsymb    = ZERO
      taux     = ZERO

      call MPI_INIT(infompi)
      comm = MPI_COMM_WORLD
      call MPI_COMM_SIZE(comm,nproc,infompi)
      call MPI_COMM_RANK(comm,me,infompi)
*
*
      tset = MPI_Wtime()
      call init_dgmres(icntl,cntl)
      if (me.eq.0) then
          call readparam_gmres(icntl,cntl,iconf,conf,restrt,Job_id)
          print *,' Parameters loaded '
          print *,' sizes   x y z',iconf(1),iconf(2),iconf(3)
          print *,' domains x y z',iconf(4),iconf(5),iconf(6)
          print *,' problem , an ',iconf(7),iconf(8)
          print *,' ac, ad, af   ',conf(3),conf(4),conf(5)
          print *,' Convection -p-',conf(1)
          print *,' Epsilon -psi- ',conf(2)
      endif
*
      call MPI_BCAST(iconf,8,MPI_INTEGER,0,comm,infompi)
      call MPI_BCAST(conf,5,MPI_DOUBLE_PRECISION,0,comm,infompi)
      call MPI_BCAST(icntl,8,MPI_INTEGER,0,comm,infompi)
      call MPI_BCAST(cntl,5,MPI_DOUBLE_PRECISION,0,comm,infompi)
      call MPI_BCAST(restrt,1,MPI_INTEGER,0,comm,infompi)
      call MPI_BCAST(Job_id,1,MPI_INTEGER,0,comm,infompi)
*
      scaling   = 0
      size_x    = iconf(1)
      size_y    = iconf(2)
      size_z    = iconf(3)
      x_domains = iconf(4)
      y_domains = iconf(5)
      z_domains = iconf(6)
      problem   = iconf(7)
      an        = iconf(8)       ! isotropic if an <> 0
      p1        = conf(1)
      p2        = conf(2)
      ac        = conf(3)
      ad        = conf(4)
      af        = conf(5)
      n         = size_x*size_y*size_z
      nxy       = (x_domains-1)*((y_domains-1)*(size_y-1)+size_y) +
     &            (y_domains-1)*((x_domains-1)*(size_x-1)+size_x) -
     &            (x_domains-1)*(y_domains-1)
      nglob     = (z_domains-1)*(x_domains*size_x - (x_domains-1))*
     &                         (y_domains*size_y - (y_domains-1)) +
     &                         z_domains*(size_z-2)*nxy + 2*nxy 
*
      ALLOCATE(abis(7*n))
      ALLOCATE(iabis(7*n))
      ALLOCATE(jabis(7*n))
      ALLOCATE(a(7*n))
      ALLOCATE(ia(7*n))
      ALLOCATE(ja(7*n))
      ALLOCATE(iwork(n))
      ALLOCATE(Perm(n))
      ALLOCATE(invPerm(n))
      ALLOCATE(x(n))
      ALLOCATE(b(n))
      ALLOCATE(temp(n))
      ALLOCATE(xsol(n))
*
      sizeIntrf = getsizeInterf(me,x_domains,y_domains,z_domains, 
     &            size_x,size_y,size_z,iwork)
*
*     Compute sum of squares of interface sizes to all neighbours
      totintrf = 2*((size_x*size_y)**2 + (size_x*size_z)**2 +
     &              (size_y*size_z)**2) +
     &         4*(size_x**2 + size_y**2 + size_z**2) + 10
      totintrfbis = 2*((size_x*size_y) + (size_x*size_z) +
     &              (size_y*size_z)) +
     &         4*(size_x + size_y + size_z) + 10
      ALLOCATE(yaux(sizeIntrf))
      ALLOCATE(weight(sizeIntrf))
      ALLOCATE(list_Intrfc(sizeIntrf))
      ALLOCATE(list_Intrfc_lo(sizeIntrf))
      ALLOCATE(Index_Intrfc(2*sizeIntrf+100))
      ALLOCATE( abb(5*sizeIntrf))
      ALLOCATE(iabb(5*sizeIntrf))
      ALLOCATE(jabb(5*sizeIntrf))
      ALLOCATE( aib(3*sizeIntrf))
      ALLOCATE(iaib(3*sizeIntrf))
      ALLOCATE(jaib(3*sizeIntrf))
      ALLOCATE( abi(3*sizeIntrf))
      ALLOCATE(iabi(3*sizeIntrf))
      ALLOCATE(jabi(3*sizeIntrf))
*
*     Allocate the dense Schur complement matrices
      ALLOCATE(   S(sizeIntrf*sizeIntrf))
      ALLOCATE(SBar(sizeIntrf*sizeIntrf))
      ALLOCATE(sqdiag_S(sizeIntrf))

*
*     Allocate the workspace for the preconjugate gradient driver
      lwork=restrt**2+(restrt*(sizeIntrf+5))+6*sizeIntrf+restrt+3
      ALLOCATE(workaux(lwork))
      ALLOCATE(orthvect(restrt))
*
*     Attach a buffer big enough for assembling the schur complement, in bytes
      max_neigintrf=max((size_x*size_y),(size_x*size_z),(size_y*size_z))
      ALLOCATE(buff(max_neigintrf**2))
      ALLOCATE(combuff(totintrf + 26*MPI_BSEND_OVERHEAD))
      sizebuff = (totintrf + 26*MPI_BSEND_OVERHEAD)*8
      call MPI_BUFFER_ATTACH(combuff, sizebuff, infompi)
*
*     Setup the interface data structure
      call setupInterf(me,x_domains,y_domains,z_domains,size_x, size_y,
     &              size_z,nbvi,indexVi,ptr_Index_Intrfc, Index_Intrfc,
     &              list_Intrfc, sizeIntrf_lo, list_Intrfc_lo,
     &              perm, invPerm, weight,iwork)
*
*     Discretise the equation (variable coefficients) on subdomains
      call local_Discret(me,x_domains,y_domains,z_domains,size_x,
     &        size_y,size_z, a, ia, ja, nnz, p1, p2, iwork,
*     &        size_y,size_z, abis, iabis, jabis, nnz, p1, p2, iwork,
     &           problem,an,ac,ad,af,Right,Left,Top,Bottom,Up,Down)

       call setup_totIntrf(me, nproc,sizeIntrf, sizeIntrf_lo,n, 
     &            totalInterf,totalInt, comm)

*      goto 555
*################################################################################
*  Put on local_Discret the a ia ja instead of abis iabis jabis/ DELETE ALLOCATE
*################################################################################
*     stupid symetrisation of the coordinate matrix a
*      k=0
*      do i=1,nnz
*        if (iabis(i).gt.jabis(i)) then
*           k=k+1
*           ia(k) = iabis(i)
*           ja(k) = jabis(i)
*           a(k)  = abis(i)
*           k=k+1
*           ia(k) = jabis(i)
*           ja(k) = iabis(i)
*           a(k)  = abis(i)
*        else if(iabis(i).eq.jabis(i)) then
*           k=k+1
*           ia(k) = iabis(i)
*           ja(k) = jabis(i)
*           a(k)  = abis(i)
*        endif
*      enddo
*########################################################
*########################################################
*########################################################
*     to test the symetry of the matrix
*########################################################
*      iabis(:)=0
*      jabis(:)=0
*      abis(:) =0.0d0
*      k=0
*      do i=1,nnz
*      do l=1,k
*         if ((ia(i).eq.jabis(l)).and.(ja(i).eq.iabis(l))) then
*             if(a(i).ne.abis(l)) print*,'ERROR on PROC',me
*             goto 333
*         endif
*      enddo
*      k=k+1
*      iabis(k) = ia(i)
*      jabis(k) = ja(i)
*      abis(k)  = a(i)
*333   enddo
*########################################################
*########################################################

!      if (me .eq. 0) print *, 'local discretisation done '
!      print *, me,' My size Interf ', sizeIntrf_lo, sizeIntrf
*********************************************************************
*         fileout=' ';
*        if(me.eq.73) then
*         write(fileout,'("MAT",I2)') me+1
*         open(unit=17,file=fileout,status="REPLACE",iostat=ios)
*         write(unit=17,FMT='(A)')"%%MatrixMarket matrix coordinate real
*     &         general" 
*         write(unit=17,FMT='(A,I4,A,I10)') "%interface size of domain ",
*     &         me+1," is:",sizeIntrf
*       write(unit=17,FMT='(I6,I6,I6)') n,n,nnz
*       write(unit=17,FMT='(I10,I10,F20.10)') (ia(i),ja(i),a(i),i=1,nnz) 
*       close(unit=17)
*        endif
*       goto 555
**********************************************************************
*         fileout=' ';
*         write(fileout,'("Invperm",I2)')me+1
*          open(unit=9,file=fileout,status="REPLACE",iostat=ios)
*        write(unit=9,FMT='(A)')"%%MatrixMarket matrix coordinate real
*     &      general" 
*        write(unit=9,FMT='(A,I4,A,I10)')"%interface size of domain ",
*     &      me+1," is:",sizeIntrf
*            write(unit=9,FMT='(I6,I6,I6)')n,1,n
*            write(unit=9,FMT='(I10,I10,I10)')(i,1,invPerm(i),i=1,n)
*            close(unit=9)
**********************************************************************

*     Permute the matrix according to the "Schur ordering"
      do i=1,nnz
         ia(i) = invPerm(ia(i))
         ja(i) = invPerm(ja(i))
      enddo
*********************************************************************
*         fileout=' ';
*         write(fileout,'("newAP",I2)') me+1
*         open(unit=18,file=fileout,status="REPLACE",iostat=ios)
*         write(unit=18,FMT='(A)')"%%MatrixMarket matrix coordinate real
*     &         general" 
*         write(unit=18,FMT='(A,I4,A,I10)') "%interface size of domain ",
*     &         me+1," is:",sizeIntrf
*       write(unit=18,FMT='(I6,I6,I6)') n,n,nnz
*       write(unit=18,FMT='(I10,I10,F15.3)') (ia(i),ja(i),a(i),i=1,nnz) 
*       close(unit=18)
*********************************************************************
*     Build the submatrices involved in the computation of the local
*     Schur complement
      call extractSubMatrices(sizeIntrf, n,
     &    a, ia, ja, nnz,
     &    abb, iabb, jabb, nnzabb, aib, iaib, jaib, nnzaib,
     &    abi, iabi, jabi, nnzabi)
*
*     Initialise the right hand side so that the solution is known
*     on the whole domain
      x(:)=1.0d0
*      x(1)=1.0d0
*      do i=2,n
*        x(i) = x(i-1)+ 1.0d0
*      enddo
*      call updateBound(x,nbvi,indexVi,
*     &      ptr_Index_Intrfc,Index_Intrfc, invPerm, buff,comm)
      call gemvcoo(n,x,b,a,ia,ja,nnz)
*
*     choose the max or the min proc of sizeIntrf to be writen on the statistique
      myIntrf(1) = sizeIntrf
      myIntrf(2) = me
      call MPI_ALLREDUCE(myIntrf,maxminIntrf(1,1),1,
     &      MPI_2INTEGER, MPI_MAXLOC,comm,infompi)
      call MPI_ALLREDUCE(myIntrf,maxminIntrf(1,2),1,
     &      MPI_2INTEGER, MPI_MINLOC,comm,infompi)
*
*      The matrix maxminIntrf(2,2) is :
*            /max_Intrf : min_Intrf/
*           /    me     :   me    /
*
      if(me.eq.maxminIntrf(2,1)) then
         write(unit=6,FMT='(A)') '@@@@@@@  sizeIntrf buffsiz_buff NBPROC
     & @@@@@@@@'
         write(unit=6,FMT='(I14,I14,I8)') sizeIntrf,totintrf,nproc
         write(unit=6,FMT='(A)') '@@@@@@@ sizeIntrf tot_Interface tot_in
     &terior maxintrf @@@@@@@'
         write(unit=6,FMT='(I14,I14,I14,I8)') sizeIntrf,totalInterf,
     &   totalInt,maxminIntrf(2,1)
       endif
     
*
      tset = MPI_Wtime() - tset

      ttot = MPI_Wtime()
*
*     Build the local Schur complement matrix
*     Define the communicator for the MUMPS package
      mumps_seq%COMM = MPI_COMM_SELF
*     Initialize an instance of the package
*     for L U factorization (sym = 0, with working host)
      mumps_seq%JOB = -1
      mumps_seq%SYM =  0
      mumps_seq%PAR =  1
      CALL DMUMPS(mumps_seq)
*     Get the statistic from the domain that have the largest interface
      if (me.ne.maxminIntrf(2,1)) mumps_seq%ICNTL(3) = 0
*      if (me.ne.0) mumps_seq%ICNTL(3) = 0
*      mumps_seq%ICNTL(3) = me+100
       mumps_seq%ICNTL(7) = 5
*      mumps_seq%ICNTL(1) = 0
*      mumps_seq%ICNTL(2) = 0
*      mumps_seq%ICNTL(3) = 0
*      mumps_seq%ICNTL(4) = 0
**      mumps_seq%write_problem = "Matrix_MUMPS_50"
      mumps_seq%N    = n
      mumps_seq%NZ   = nnz
      mumps_seq%IRN  => ia
      mumps_seq%JCN  => ja
      mumps_seq%A    => a
*      mumps_seq%RHS  => b
      mumps_seq%SIZE_SCHUR = sizeIntrf
      if ( sizeIntrf.ne.0) then
         ALLOCATE( mumps_seq%LISTVAR_SCHUR(sizeIntrf) )
         mumps_seq%SCHUR => S
         do i=1,sizeIntrf 
           mumps_seq%LISTVAR_SCHUR(i) = i
         enddo
         mumps_seq%ICNTL(19) = 1
      else
         mumps_seq%ICNTL(19) = 0
      endif
*
*     Build the local Schur complement using MUMPS
*     symbolic analysis
      mumps_seq%JOB       = 1
      tsymb = MPI_Wtime() !call secdeb(tsymb)
      CALL DMUMPS(mumps_seq)
      tsymb = MPI_Wtime()-tsymb  ! call secfin(tsymb)

*     Numerical factorization
      call mpi_barrier(comm,infompi)
      mumps_seq%JOB       = 2
      tfacto = MPI_Wtime()!call secdeb(tfacto)
      CALL DMUMPS(mumps_seq)
*      print *, me, 'Numerical Factorization done'
      if (mumps_seq%INFO(1).ne.0) then
        print *,'***** Error on proc ', me, mumps_seq%INFO(1)
        goto 444
        mumps_seq%JOB = -2
        CALL DMUMPS(mumps_seq)
        call MPI_ABORT(infompi)
        call MPI_FINALIZE(infompi)
        stop
      endif
      tfacto = MPI_Wtime()-tfacto !  call secfin(tfacto)
*
*     Setup the right hand side of the Schur system
      do i=1,sizeIntrf
        temp(i) = ZERO
      enddo
*     temp = bI
      do i=sizeIntrf+1,n
         temp(i) = b(i)
      enddo
*     temp = AII^{-1}*bI
      mumps_seq%RHS  => temp
      mumps_seq%JOB       = 3
      tsolv = MPI_Wtime() !call secdeb(tsolv)
      CALL DMUMPS(mumps_seq)
      tsolv = MPI_Wtime()- tsolv !call secfin(tsolv)
*************************************************************************************
* faut faire attention a gemvcoo car la multiplication se fait suivant ia(i)
* et ja(i) et donc suivant les indice globale de chaque matrice a du subdomain
* et donc si dans notre cas on a mis les interfaces a la fin de la matrice on aura 
* un probleme car le vecteur yaux est de taille sizeIntrf donc le mieux est d'allouer
* le vecteur resultat de taille n car dans la multiplication si les indices intefaces
* sont a la fin alors c'est yaux(n-sizeIntrf+1:n) qui sera rempli du resultat
      call mpi_barrier(comm,infompi)
      if(me.eq.0) print*,me,'===============FINISH SOLVE============='
*     yaux = ABI * AII^{-1}* bI
      call gemvcoo(sizeIntrf,temp,yaux,
     &            abi,iabi,jabi,nnzabi) ! faut faire attention a gemvcoo 
      do i=1,sizeIntrf
         b(i) = b(i) - yaux(i)
      enddo
*
*     Assemble the right hand side on the boundaries
      call updateBound(b,nbvi,indexVi,
     &      ptr_Index_Intrfc,Index_Intrfc, invPerm, buff,comm)
      if(me.eq.0) print*,me,'=============FINISH UPDATE============='
*
*     Initialisation for the coarse grid 
*      nbdomain = x_domains*y_domains*z_domains
*      ALLOCATE(Acoa(nbdomain,nbdomain))
*      ALLOCATE(Rt(sizeintrf,nbvi+1))
*      ALLOCATE(tempAcR(nbdomain))
*      ALLOCATE(tempAcRbis(nbdomain))
*      ALLOCATE(IPIV(nbdomain))
*      Acoa(:,:)=0.0d0
*      Rt(:,:)=0.0d0
**     Added on the last value of indexVi to use for the coarse grid
*      indexVi(nbvi+1)=me
**
*       if(me.eq.0)print*,me,'=============START COARSE============='
*     
*      call CoarsePcondprep(Acoa,Rt,S,sizeIntrf,nbdomain,nbvi,weight,
*     &      indexVi,ptr_Index_Intrfc,Index_Intrfc,invPerm,IPIV,comm)
*
*      if(me.eq.0)print*,me,'==========FINISH COARSE ============'

*     used for the estimation of the condition number
      ALLOCATE(wor(4*sizeIntrf))
      ALLOCATE(iwor(sizeIntrf))
*     Copy and then assemble the local Schur complement
      call mpi_barrier(comm,infompi)
*     ###############################################
*       ATTENTION MUMPS give the Transpose of Schur
*     ###############################################
      tpreloc = MPI_Wtime() !call secdeb(tpreloc)
      tassemb = MPI_Wtime() 
      call DCOPY(sizeIntrf*sizeIntrf,S,1,SBar,1)
      call assembleSchur(SBar,sizeIntrf,nbvi,indexVi,
     &                   ptr_Index_Intrfc, Index_Intrfc,
     &                   invPerm,buff)

      tassemb = MPI_Wtime() -tassemb
      call mpi_barrier(comm,infompi)
      if(me.eq.0)print*,me,'==========FINISH ASSEMBLE ============'
*     scale the Schur complement and the precoditionner 
*     first compute the scaling matrix that is the diagonal 
*     of the Schur complement(use the assembled Schur)
*     then scale the local Schur matrix (assembled and non assembled)
      if(scaling.eq.1) then
         call scale_S(SBar,SizeIntrf,sqdiag_S,'N')
         call scale_S(S,SizeIntrf,sqdiag_S,'Y')
      endif
*     compute the "xx" norm of the Schur complement 
C      Normm= DLANGE( '1', sizeIntrf,sizeIntrf,SBar,sizeIntrf,WOR)
      call mpi_barrier(comm,infompi)
      tprefact = MPI_Wtime() 
*     Prefactorize (L U) the locally assembled Schur complement
      ALLOCATE(IPIV_S(sizeIntrf))
      call DGETRF( sizeIntrf,sizeIntrf,SBar,sizeIntrf,IPIV_S,ierr)
      if(ierr.ne.0) then
         print *,me, '########################################'
         print *,me, '           ERROR Pcond Facto  ',ierr
         print *,me, '########################################'
         goto 444
      endif
      tprefact = MPI_Wtime()-tprefact
      tpreloc = MPI_Wtime()-tpreloc !call secfin(tpreloc)
      if(me.eq.0)print*,me,'=======FINISH Facto Pcond =========='
*     Estimate the condition number of the assembled Schur complement in "Normm" norm
C      call DGECON('1',sizeIntrf,SBar,sizeIntrf,Normm,
C     &                           Rcond,wor,iwor,ierr)
*
*     scale the RHS "b" by the diagonal of the  Schur complement
*     i.e b =Diag{-1/2}*b in the case of scaling=1
      if(scaling.eq.1) then
         do i=1,sizeIntrf
            b(i)=b(i)*sqdiag_S(i)
         enddo
      endif
*
       if(me.eq.0)print*,me,'========= START GMRES ============='
*
*     Initialize and start the Gradient conjugate solver
      do i=1,sizeIntrf
         workaux(sizeIntrf+i) = b(i)
      enddo
*     Setup the initial guess
*      do i=1,sizeIntrf
*        workaux(i) = 1.0d0
*      enddo
      if (me.eq.maxminIntrf(2,1)) then
        icntl(1) = 6
        icntl(2) = 6
        icntl(3) = 6
      else
        icntl(1) = 0
        icntl(2) = 0
        icntl(3) = 0
      endif
      call mpi_barrier(comm,infompi)
      tcg = MPI_Wtime() !call secdeb(tcg)
10    call drive_dgmres(totalInterf,sizeIntrf,restrt,lwork,workaux,irc,
     &             icntl,cntl,info,rinfo)
       revcom = irc(1)
       colx   = irc(2)
       coly   = irc(3)
       colz   = irc(4)
       nbscal = irc(5)
*
*     ###############################################
*       ATTENTION MUMPS give the Transpose of Schur
*     ###############################################
       if (revcom.eq.1) then
*         perform the matrix vector product
*         workaux(colz) <-- S * workaux(colx)
*          if (me.eq.maxminIntrf(2,1))
*     &                         print*,'perform Matrix vector product'
          taux = MPI_Wtime() 
          call DGEMV('T',sizeIntrf,sizeIntrf,ONE,S,
     &       sizeIntrf,workaux(colx),1, ZERO, workaux(colz),1)
          tMatVect = MPI_Wtime()-taux+tMatVect
          call updateBound(workaux(colz),nbvi,indexVi,
     &       ptr_Index_Intrfc,Index_Intrfc,invPerm, buff,comm)
C          tMatVect = MPI_Wtime()-taux+tMatVect
C          print*,me,'------------- fin MATVECT product ------------'
          goto 10
       else if (revcom.eq.2) then
          taux = MPI_Wtime() 
          call pcondLocSnosym(workaux(colx),sizeIntrf,
     &         workaux(colz),SBar,IPIV_S,ierr)  
          if(ierr.ne.0) goto 556
C          tpreapp = MPI_Wtime()-taux+tpreapp
          call updateBound(workaux(colz),nbvi,indexVi,
     &         ptr_Index_Intrfc,Index_Intrfc,invPerm, buff,comm)
          tpreapp = MPI_Wtime()-taux+tpreapp
          goto 10
       else if (revcom.eq.3) then
*         perform the right preconditionning workaux(colz) <-- M^{-1} * workaux(colx)
*          if (me.eq.maxminIntrf(2,1)) print*,'perform precond solve'
          taux = MPI_Wtime() 
          call pcondLocSnosym(workaux(colx),sizeIntrf,
     &         workaux(colz),SBar,IPIV_S,ierr)  
          if(ierr.ne.0) goto 556
C          tpreapp = MPI_Wtime()-taux+tpreapp
          call updateBound(workaux(colz),nbvi,indexVi,
     &         ptr_Index_Intrfc,Index_Intrfc,invPerm, buff,comm)
          tpreapp = MPI_Wtime()-taux+tpreapp
*          perform the coarse preconditionner
*          call CoarsePcondsolv(Acoa,Rt,workaux(colx),workaux(colz),yaux,
*     &     tempAcR,tempAcRbis,sizeIntrf,nbdomain,nbvi,indexVi,IPIV,comm)
C          tpreapp = MPI_Wtime()-taux+tpreapp
          goto 10
       else if (revcom.eq.4) then
*         perform the scalar product workaux(colz) <-- workaux(colx)^T * weight .* workaux(coly)
*          if (me.eq.maxminIntrf(2,1)) print*,'perform dot product'
          taux = MPI_Wtime() 
          call weighted(sizeIntrf,workaux(coly),weight,yaux) 
          call DGEMV('T',sizeIntrf,nbscal,ONE,workaux(colx),
     &         sizeIntrf,yaux,1,ZERO,orthvect,1)
          call MPI_ALLREDUCE(orthvect,workaux(colz),nbscal,
     &      MPI_DOUBLE_PRECISION, MPI_SUM,comm,infompi)
          tdot = MPI_Wtime()-taux+tdot
          goto 10
       endif
      tcg = MPI_Wtime()-tcg !call secfin(tcg)
*
       res(1) = ddotwS(sizeIntrf,list_Intrfc,
     &          workaux(1), workaux(1),weight,invPerm) 
       call MPI_ALLREDUCE(res,res(2),1,
     &      MPI_DOUBLE_PRECISION, MPI_SUM,comm,infompi)
*
*     return to the reel solution "u", i.e u =Diag{-1/2}*y 
*     where y is the solution of the GM in the case of scaling=1
      if(scaling.eq.1) then
         do i=1,sizeIntrf
            workaux(i) = workaux(i) * sqdiag_S(i)
         enddo
      endif
*
*      End of the GM solution phase
*     ****************************************************************
*                perform the interior solution of the system
*     ****************************************************************
      do i=sizeIntrf+1,n
        temp(i) = ZERO
      enddo
*     temp = xB
      do i=1,sizeIntrf
         temp(i) = workaux(i)
      enddo
      do i=1,sizeIntrf
       xsol(i) = ZERO
      enddo
*     xaux = AIB * xB
      call gemvcoo(n-sizeIntrf,temp,xsol,
     &            aib,iaib,jaib,nnzaib)
*     Note that xsol(n-sizeIntrf+1:n) contient le resultat de cette multiplication
*     c'est pour cela apres je met xsol(i)
*     temp=bI-AIB*xB
      do i=sizeIntrf+1,n
         temp(i) = b(i) - xsol(i)
      enddo
      do i=1,sizeIntrf
       temp(i) = ZERO
      enddo
*     temp=AII{-1}*(bI-AIB*xB)
      mumps_seq%RHS  => temp
      mumps_seq%JOB       = 3
      tsolvi = MPI_Wtime()!call secdeb(taux)
      CALL DMUMPS(mumps_seq)
      tsolvi = MPI_Wtime()-tsolvi !call secfin(taux)
      tsolv = tsolv + tsolvi
*     perform the solution xsol on the local subdomain
      do i=1,sizeIntrf
        xsol(i) = workaux(i)
      enddo
      do i=sizeIntrf+1,n
         xsol(i) = temp(i)
      enddo
      ttot = MPI_Wtime()-ttot

**********************************************************************
* write the local solution of all entries on a file
**********************************************************************
*      if (me.eq.0) then
*         fileout='';
*         write(fileout,'("XSOLUTION",I4)'),me+1
*          open(unit=31,file="XSOL",status="REPLACE",iostat=ios)
*          if (ios.ne.0) then
*             print *,"ATTENTION unable to open file XSOLUTION"
*             goto 555
*          else 
*            write(unit=31,'(A)'),"%%MatrixMarket matrix coordinate real
*     &      general" 
*            write(unit=31,'(A,I4,A,I10)'),"%interface size of domain ",
*     &      me+1," is:",sizeIntrf
*            write(unit=31,'(I6,I6,I6)'),n,1,n
*            write(unit=31,'(I10,I10,F15.3)'),(i,1,xsol(i),i=1,n)
*            close(unit=31)
*          endif
*      endif 

**********************************************************************
      if (me.eq.maxminIntrf(2,1)) then
        do i = 1, min(50, n)
          print *,me, i, xsol(i)
        enddo
      endif
**********************************************************************
*     ***************************
*            Norm of Solution
*     ***************************
*     interface norm
      Norm_intrf=0.0d0
      Norm_int  =0.0d0
      DO i=1,sizeIntrf
         Norm_intrf = Norm_intrf + weight(i)*(xsol(i)**2)
      ENDDO
*     interior norm
      Norm_int =  DDOT(n-sizeIntrf,xsol(sizeIntrf+1),1,
     &                 xsol(sizeIntrf+1),1)
*     all vector reduce sum and sqrt    
      Normmsol    =  Norm_intrf + Norm_int
      CALL MPI_ALLREDUCE(Normmsol,Normmbis, 1,  
     &     MPI_DOUBLE_PRECISION, MPI_SUM,comm,infompi)
      Normmsol = SQRT(Normmbis)
      if(me.eq.maxminIntrf(2,1))write(unit=6,FMT='(A30,1P,E26.18)')
     & '=====> Voici norm solution =====> ',Normmsol
*     ***************************
**********************************************************************
*
*
*     Allocate the timings arrays
*
*
      ALLOCATE(o_tot(nproc))
      ALLOCATE(o_tot_bis(nproc))
      ALLOCATE(o_tsymb(nproc))
      ALLOCATE(o_tsymb_bis(nproc))
      ALLOCATE(o_tfacto(nproc))
      ALLOCATE(o_tfacto_bis(nproc))
      ALLOCATE(o_tsolv(nproc))
      ALLOCATE(o_tsolv_bis(nproc))

      ALLOCATE(o_tpreloc(nproc))
      ALLOCATE(o_tpreloc_bis(nproc))
      ALLOCATE(o_tassemb(nproc))
      ALLOCATE(o_tassemb_bis(nproc))
      ALLOCATE(o_tprefact(nproc))
      ALLOCATE(o_tprefact_bis(nproc))

      ALLOCATE(o_tcg(nproc))
      ALLOCATE(o_tcg_bis(nproc))
      ALLOCATE(o_tMatVect(nproc))
      ALLOCATE(o_tMatVect_bis(nproc))
      ALLOCATE(o_tpreapp(nproc))
      ALLOCATE(o_tpreapp_bis(nproc))
      ALLOCATE(o_tdot(nproc))
      ALLOCATE(o_tdot_bis(nproc))

      ALLOCATE(o_memo(nproc))
      ALLOCATE(o_memo_bis(nproc))
      ALLOCATE(o_effmemo(nproc))
      ALLOCATE(o_effmemo_bis(nproc))
      ALLOCATE(o_schur(nproc))
      ALLOCATE(o_schur_bis(nproc))

      do i=1,nproc
        o_tot(i)         = ZERO
        o_tot_bis(i)     = ZERO
        o_tsymb(i)       = ZERO
        o_tsymb_bis(i)   = ZERO
        o_tfacto(i)      = ZERO
        o_tfacto_bis(i)  = ZERO
        o_tsolv(i)       = ZERO
        o_tsolv_bis(i)   = ZERO

        o_tpreloc(i)     = ZERO
        o_tpreloc_bis(i) = ZERO
        o_tassemb        = ZERO
        o_tassemb_bis    = ZERO
        o_tprefact       = ZERO
        o_tprefact_bis   = ZERO

        o_tcg(i)         = ZERO
        o_tcg_bis(i)     = ZERO
        o_tMatVect(i)    = ZERO
        o_tMatVect_bis(i)= ZERO
        o_tpreapp(i)     = ZERO
        o_tpreapp_bis(i) = ZERO
        o_tdot(i)        = ZERO
        o_tdot_bis(i)    = ZERO

        o_memo(i)        = 0
        o_memo_bis(i)    = 0
        o_effmemo(i)     = 0
        o_effmemo_bis(i) = 0
        o_schur(i)       = 0
        o_schur_bis(i)   = 0

      enddo
*
*     Get Timings array

      o_tsymb(me+1)    = tsymb
      o_tfacto(me+1)   = tfacto
      o_tsolv(me+1)    = tsolv
      o_tpreloc(me+1)  = tpreloc
      o_tassemb(me+1)  = tassemb
      o_tprefact(me+1) = tprefact
      o_tcg(me+1)      = tcg
      o_tMatVect(me+1) = tMatVect
      o_tpreapp(me+1)  = tpreapp
      o_tdot(me+1)     = tdot
      o_schur(me+1)    = sizeIntrf
      o_memo(me+1)     = mumps_seq%INFOG(19)
      o_effmemo(me+1)  = mumps_seq%INFOG(22)
      o_tot(me+1)      = ttot
*
      call MPI_REDUCE(o_tot,o_tot_bis,nproc,
     &      MPI_DOUBLE_PRECISION, MPI_SUM,maxminIntrf(2,1),comm,infompi)
      call MPI_REDUCE(o_tsymb,o_tsymb_bis,nproc,
     &      MPI_DOUBLE_PRECISION, MPI_SUM,maxminIntrf(2,1),comm,infompi)
      call MPI_REDUCE(o_tfacto,o_tfacto_bis,nproc,
     &      MPI_DOUBLE_PRECISION, MPI_SUM,maxminIntrf(2,1),comm,infompi)
      call MPI_REDUCE(o_tsolv,o_tsolv_bis,nproc,
     &      MPI_DOUBLE_PRECISION, MPI_SUM,maxminIntrf(2,1),comm,infompi)
      call MPI_REDUCE(o_tpreloc,o_tpreloc_bis,nproc,
     &      MPI_DOUBLE_PRECISION, MPI_SUM,maxminIntrf(2,1),comm,infompi)
      call MPI_REDUCE(o_tassemb,o_tassemb_bis,nproc,
     &      MPI_DOUBLE_PRECISION, MPI_SUM,maxminIntrf(2,1),comm,infompi)
      call MPI_REDUCE(o_tprefact,o_tprefact_bis,nproc,
     &      MPI_DOUBLE_PRECISION, MPI_SUM,maxminIntrf(2,1),comm,infompi)
      call MPI_REDUCE(o_tcg,o_tcg_bis,nproc,
     &      MPI_DOUBLE_PRECISION, MPI_SUM,maxminIntrf(2,1),comm,infompi)
      call MPI_REDUCE(o_tMatVect,o_tMatVect_bis,nproc,
     &      MPI_DOUBLE_PRECISION, MPI_SUM,maxminIntrf(2,1),comm,infompi)
      call MPI_REDUCE(o_tpreapp,o_tpreapp_bis,nproc,
     &      MPI_DOUBLE_PRECISION, MPI_SUM,maxminIntrf(2,1),comm,infompi)
      call MPI_REDUCE(o_tdot,o_tdot_bis,nproc,
     &      MPI_DOUBLE_PRECISION, MPI_SUM,maxminIntrf(2,1),comm,infompi)
      call MPI_REDUCE(o_schur,o_schur_bis,nproc,
     &      MPI_INTEGER, MPI_SUM,maxminIntrf(2,1),comm,infompi)
      call MPI_REDUCE(o_memo,o_memo_bis,nproc,
     &      MPI_INTEGER, MPI_SUM,maxminIntrf(2,1),comm,infompi)
      call MPI_REDUCE(o_effmemo,o_effmemo_bis,nproc,
     &      MPI_INTEGER, MPI_SUM,maxminIntrf(2,1),comm,infompi)
*
*     collecte the small and the big condition number for statistic
      call MPI_ALLREDUCE(Rcond,Rcondmax,1,MPI_DOUBLE_PRECISION,
     &                                    MPI_MAX,comm,infompi)
      call MPI_ALLREDUCE(Rcond,Rcondmin,1,MPI_DOUBLE_PRECISION,
     &                                    MPI_MIN,comm,infompi)
*
      avertot    = ZERO
      aversymb   = ZERO
      averfact   = ZERO
      aversolv   = ZERO
      averpreloc = ZERO
      averassemb = ZERO
      averprefact= ZERO
      avercg     = ZERO
      averMatVect= ZERO
      averpreapp = ZERO
      averdot    = ZERO
      avertime   = ZERO
*
      mintot    = 9999999d0
      minsymb   = 9999999d0
      minfact   = 9999999d0
      minsolv   = 9999999d0
      minpreloc = 9999999d0
      minassemb = 9999999d0
      minprefact= 9999999d0
      mincg     = 9999999d0
      minMatVect= 9999999d0
      minpreapp = 9999999d0
      mindot    = 9999999d0
      mintime   = 9999999d0
*
      maxtot    = ZERO
      maxsymb   = ZERO
      maxfact   = ZERO
      maxsolv   = ZERO
      maxpreloc = ZERO
      maxassemb = ZERO
      maxprefact= ZERO
      maxcg     = ZERO
      maxMatVect= ZERO
      maxpreapp = ZERO
      maxdot    = ZERO
      maxtime   = ZERO
*
      minmemo      = 9999999
      maxmemo      = 0
      avermemo     = 0
      mineffmemo   = 9999999
      maxeffmemo   = 0
      avereffmemo  = 0
      minschur     = 9999999
      maxschur     = 0
      averschur    = 0

      if (me .eq.maxminIntrf(2,1)) then
        do i = 1, nproc
          avertot    = avertot + o_tot_bis(i)
          mintot    = min(mintot , o_tot_bis(i) )
          maxtot    = max(maxtot , o_tot_bis(i) )

          aversymb   = aversymb + o_tsymb_bis(i)
          averfact   = averfact + o_tfacto_bis(i)
          aversolv   = aversolv + o_tsolv_bis(i)
          averpreloc = averpreloc + o_tpreloc_bis(i)
          averassemb = averassemb + o_tassemb_bis(i)
          averprefact= averprefact + o_tprefact_bis(i)
          avercg     = avercg + o_tcg_bis(i)
          averMatVect= averMatVect + o_tMatVect_bis(i)
          averpreapp = averpreapp + o_tpreapp_bis(i)
          averdot    = averdot + o_tdot_bis(i)

          minsymb   = min(minsymb , o_tsymb_bis(i) )
          minfact   = min(minfact , o_tfacto_bis(i) )
          minsolv   = min(minsolv , o_tsolv_bis(i) )
          minpreloc = min(minpreloc , o_tpreloc_bis(i) )
          minassemb = min(minassemb , o_tassemb_bis(i) )
          minprefact= min(minprefact , o_tprefact_bis(i) )
          mincg     = min(mincg , o_tcg_bis(i) )
          minMatVect= min(minMatVect , o_tMatVect_bis(i) )
          minpreapp = min(minpreapp , o_tpreapp_bis(i) )
          mindot    = min(mindot , o_tdot_bis(i) )

          maxsymb   = max(maxsymb , o_tsymb_bis(i) )
          maxfact   = max(maxfact , o_tfacto_bis(i) )
          maxsolv   = max(maxsolv , o_tsolv_bis(i) )
          maxpreloc = max(maxpreloc , o_tpreloc_bis(i) )
          maxassemb = max(maxassemb , o_tassemb_bis(i) )
          maxprefact= max(maxprefact , o_tprefact_bis(i) )
          maxcg     = max(maxcg , o_tcg_bis(i) )
          maxMatVect= max(maxMatVect , o_tMatVect_bis(i) )
          maxpreapp = max(maxpreapp , o_tpreapp_bis(i) )
          maxdot    = max(maxdot , o_tdot_bis(i) )

          avermemo      = avermemo + o_memo_bis(i)
          avereffmemo   = avereffmemo + o_effmemo_bis(i)
          minmemo       = min(minmemo , o_memo_bis(i) )
          maxmemo       = max(maxmemo , o_memo_bis(i) )
          mineffmemo    = min(mineffmemo , o_effmemo_bis(i) )
          maxeffmemo    = max(maxeffmemo , o_effmemo_bis(i) )
          averschur     = averschur + o_schur_bis(i)
          minschur      = min(minschur , o_schur_bis(i) )
          maxschur      = max(maxschur , o_schur_bis(i) )

        enddo
        avermemo   = int(avermemo/nproc)
        avereffmemo= int(avereffmemo/nproc)
        averschur   = int(averschur/nproc)

        avertot    = avertot/nproc 
        aversymb   = aversymb/nproc
        averfact   = averfact/nproc
        aversolv   = aversolv/nproc
        averpreloc = averpreloc/nproc
        averassemb = averassemb/nproc
        averprefact= averprefact/nproc
        avercg     = avercg/nproc 
        averMatVect= averMatVect/nproc
        averpreapp = averpreapp/nproc 
        averdot    = averdot/nproc 
        avertime   = aversymb+averfact+averpreloc+avercg+aversolv   
        ttotal     = tcg+tpreloc+tfacto+tsymb+tsolv
*
        write(*,*)
        write(*,*)
        write(*,*)
        if ((icntl(7).eq.1).and.(info(2).gt.2)) then
          write(*,*)
          write(*,'(A19,E25.16E3)')  '  M_AS Approx. Cond :     ',
     &    rinfo(2)/rinfo(3)
          write(*,'(A19,E25.16E3,A5,E25.16E3)') 
     &    '  M_AS Eig min Max  :     ',rinfo(3),'---',rinfo(2)
        endif
C        write(unit=6,FMT='(A,1P,E26.18,A9,1P,E26.18)')
C     &        ' M_AS The m/M Pcond Condition number   ',
C     &        1.0/Rcondmax,'--Azz--',1.0/Rcondmin

        write(*,*)
        write(unit=6,FMT='(A,I10)')
     &        ' M_AS The number of Processors             = ',nproc
        if (info(1).eq.-3)then
        write(unit=6,FMT='(A,A10)')
     &        ' M_AS The number of iterations M_AS        = ','000'
        else
        write(unit=6,FMT='(A,I10)')
     &        ' M_AS The number of iterations M_AS        = ',info(2)
        endif

        write(unit=6,FMT='(A,I10)') 
     &  ' M_AS  total Intreface (iterative system)  = ',totalInterf

        write(unit=6,FMT='(A112)')
     &  '--------------------------------------------------------------'
        write(unit=6,FMT='(A112)')
     &  '--------------------------------------------------------------'
        write(unit=6,FMT='(A47,A10,A9,A10,A9,A10,A9,I8)') 
     &  '','MIN  ','','AVG  ','','MAX  ','',maxminIntrf(2,1)
        write(unit=6,FMT='(A112)')
     &  '--------------------------------------------------------------'

        write(unit=6,FMT='(A,I10,A9,I10,A9,I10,A9,I10)') 
     &        ' M_AS Schur size                           = ',
     &        minschur,'--Azz--',averschur,'--Azz--',
     &        maxschur,'--Azz--',o_schur(me+1)
        write(unit=6,FMT='(A,I10,A9,I10,A9,I10,A9,I10)') 
     &        ' M_AS Schur           Memory (MB)          = ',
     &        int(((minschur**2)*8)/1.0d6),'--Azz--',
     &        int(((averschur**2)*8)/1.0d6),'--Azz--',
     &        int(((maxschur**2)*8)/1.0d6),'--Azz--',
     &        int(((o_schur(me+1)**2)*8)/1.0d6)

        write(unit=6,FMT='(A,I10,A9,I10,A9,I10,A9,I10)') 
     &        ' M_AS MUMPS           Memory INFOG(19)     = ',
     &        minmemo,'--Azz--',avermemo,'--Azz--',
     &        maxmemo,'--Azz--',o_memo(me+1)
        write(unit=6,FMT='(A,I10,A9,I10,A9,I10,A9,I10)') 
     &        ' M_AS MUMPS effective Memory INFOG(22)     = ',
     &        mineffmemo,'--Azz--',avereffmemo,'--Azz--',
     &        maxeffmemo,'--Azz--',o_effmemo(me+1)
        write(unit=6,FMT='(A112)')
     &  '--------------------------------------------------------------'
        write(unit=6,FMT='(A,F10.3,A9,F10.3,A9,F10.3,A9,F10.3)') 
     &        ' M_AS MUMPS  Symbolique       time         = ',
     &        minsymb,'--Azz--',aversymb,'--Azz--',
     &        maxsymb,'--Azz--',tsymb
        write(unit=6,FMT='(A,F10.3,A9,F10.3,A9,F10.3,A9,F10.3)') 
     &        ' M_AS MUMPS  Factorization    time         = ',
     &        minfact,'--Azz--',averfact,'--Azz--',
     &        maxfact,'--Azz--',tfacto
        write(unit=6,FMT='(A,F10.3,A9,F10.3,A9,F10.3,A9,F10.3)') 
     &        ' M_AS MUMPS  Solve            time         = ',
     &        minsolv,'--Azz--',aversolv,'--Azz--',
     &        maxsolv,'--Azz--',tsolv
        write(unit=6,FMT='(A,F10.3,A9,F10.3,A9,F10.3,A9,F10.3)')
     &        ' M_AS Assemble                time         = ',
     &        minassemb,'--Azz--',averassemb,'--Azz--',
     &        maxassemb,'--Azz--',tassemb
        write(unit=6,FMT='(A,F10.3,A9,F10.3,A9,F10.3,A9,F10.3)')
     &        ' M_AS Precond Facto           time         = ',
     &        minprefact,'--Azz--',averprefact,'--Azz--',
     &        maxprefact,'--Azz--',tprefact
        write(unit=6,FMT='(A,F10.3,A9,F10.3,A9,F10.3,A9,F10.3)')
     &        ' M_AS Preloc                  time         = ',
     &        minpreloc,'--Azz--',averpreloc,'--Azz--',
     &        maxpreloc,'--Azz--',tpreloc
        write(unit=6,FMT='(A,F10.3,A9,F10.3,A9,F10.3,A9,F10.3)')
     &        ' M_AS GMRES                   time         = ',
     &        mincg,'--Azz--',avercg,'--Azz--',
     &        maxcg,'--Azz--',tcg
        write(unit=6,FMT='(A,F10.3,A9,F10.3,A9,F10.3,A9,F10.3)')
     &        ' M_AS GMRES                   time per iter= ',
     &        mincg/info(2),'--Azz--',avercg/info(2),'--Azz--',
     &        maxcg/info(2),'--Azz--',tcg/info(2)
        write(unit=6,FMT='(A,F10.3,A9,F10.3,A9,F10.3,A9,F10.3)')
     &        ' M_AS GMRES-MatVect           time         = ',
     &        minMatVect,'--Azz--',averMatVect,'--Azz--',
     &        maxMatVect,'--Azz--',tMatVect
        write(unit=6,FMT='(A,F10.3,A9,F10.3,A9,F10.3,A9,F10.3)')
     &        ' M_AS GMRES-MatVect           time per iter= ',
     &        minMatVect/(info(2)+1),'--Azz--',
     &        averMatVect/(info(2)+1),'--Azz--',
     &        maxMatVect/(info(2)+1),'--Azz--',
     &        tMatVect/(info(2)+1)
        write(unit=6,FMT='(A,F10.3,A9,F10.3,A9,F10.3,A9,F10.3)')
     &        ' M_AS GMRES-Preapply          time         = ',
     &        minpreapp,'--Azz--',averpreapp,'--Azz--',
     &        maxpreapp,'--Azz--',tpreapp
        write(unit=6,FMT='(A,F10.3,A9,F10.3,A9,F10.3,A9,F10.3)')
     &        ' M_AS GMRES-Preapply          time per iter= ',
     &        minpreapp/(info(2)+1),'--Azz--',
     &        averpreapp/(info(2)+1),'--Azz--',
     &        maxpreapp/(info(2)+1),'--Azz--',
     &        tpreapp/(info(2)+1)
        write(unit=6,FMT='(A,F10.3,A9,F10.3,A9,F10.3,A9,F10.3)')
     &        ' M_AS GMRES-dot               time         = ',
     &        mindot,'--Azz--',averdot,'--Azz--',
     &        maxdot,'--Azz--',tdot
        write(unit=6,FMT='(A112)')
     &  '--------------------------------------------------------------'
        write(unit=6,FMT='(A,F10.3,A9,F10.3,A9,F10.3,A9,F10.3)')
     &        ' M_AS  Total SUM of all time               = ',
     &        minsymb+minfact+minsolv+minassemb+minprefact+
     &        mincg,'--Azz--',
     &        aversymb+averfact+aversolv+averassemb+averprefact+
     &        avercg,'--Azz--',
     &        maxsymb+maxfact+maxsolv+maxassemb+maxprefact+
     &        maxcg,'--Azz--',
     &        tsymb+tfacto+tsolv+tassemb+tprefact+tcg
        write(unit=6,FMT='(A,F10.3,A9,F10.3,A9,F10.3,A9,F10.3)')
     &        ' M_AS   AZZAM Global computing time        = ',
     &        mintot,'--Azz--',avertot,'--Azz--',
     &        maxtot,'--Azz--',ttot
        write(unit=6,FMT='(A112)')
     &  '--------------------------------------------------------------'

        write(unit=6,FMT='(A,I10,A9,I10)')
     &        ' M_AS AZZAM MUMPS info(7) and info(8)      : ',
     &        mumps_seq%INFO(7),'--->',mumps_seq%INFO(8)
        write(unit=6,FMT='(A,F10.3)')
     &        ' M_AS AZZAM SUM  of local compute time egal: ',ttotal
        write(unit=6,FMT='(A,F10.3)')
     &        ' M_AS AZZAM GLOBAL compute time egal       : ',ttot
        write(unit=6,FMT='(A,F10.3)')
     &        ' M_AS AZZAM SET UP time egal               : ',tset
      endif
*
*        call MPI_BCAST(aversymb,1,MPI_DOUBLE_PRECISION,0,comm,infompi)
*        call MPI_BCAST(averfact,1,MPI_DOUBLE_PRECISION,0,comm,infompi)
*        call MPI_BCAST(averpreloc,1,MPI_DOUBLE_PRECISION,0,comm,infompi)
*        call MPI_BCAST(avercg,1,MPI_DOUBLE_PRECISION,0,comm,infompi)
*        call MPI_BCAST(aversolv,1,MPI_DOUBLE_PRECISION,0,comm,infompi)
*        call MPI_BCAST(avertime,1,MPI_DOUBLE_PRECISION,0,comm,infompi)

*     open a statistic file and write the statistic
      if (me.eq.maxminIntrf(2,1)) then 
       open(unit=1111,file="STAT_conv",status="OLD",position="APPEND",
     & iostat=ios)
       write(unit=1111,FMT='(I6)',ADVANCE="NO") Job_id
       write(unit=1111,FMT='(I3,A,I2,A,I2,A,I2,A,I2,A,I2,A)',
     &                                          ADVANCE="NO") 
     & x_domains,"(",size_x,")", y_domains,"(",size_y,")",
     & z_domains,"(",size_z,")"
       write(unit=1111,FMT='(I7,A1,I5)',ADVANCE="NO") maxminIntrf(1,2),
     &                                            "/",maxminIntrf(1,1)
       write(unit=1111,FMT='(1PD12.3)',ADVANCE="NO") mumps_seq%RINFO(3) 

       write(unit=1111,FMT='(F9.3)',ADVANCE="NO") 
     &          (((DBLE(mumps_seq%INFO(9))*8.0d0)+ 
     &          DBLE(mumps_seq%INFO(10))*4.0d0)/DBLE(1000000))
       write(unit=1111,FMT='(I8)',ADVANCE="NO") mumps_seq%INFO(16)

       write(unit=1111,FMT='(I6)',ADVANCE="NO") info(2)
       write(unit=1111,FMT='(F9.3)',ADVANCE="NO") tsymb
       write(unit=1111,FMT='(F9.3)',ADVANCE="NO") tfacto
       write(unit=1111,FMT='(F9.3)',ADVANCE="NO") tsolv

       write(unit=1111,FMT='(F9.3)',ADVANCE="NO") tpreloc
       write(unit=1111,FMT='(F9.3)',ADVANCE="NO") tassemb
       write(unit=1111,FMT='(F9.3)',ADVANCE="NO") tprefact

       write(unit=1111,FMT='(F9.3)',ADVANCE="NO") tcg
       write(unit=1111,FMT='(F9.3)',ADVANCE="NO") tMatVect
       write(unit=1111,FMT='(F9.3)',ADVANCE="NO") tpreapp
       write(unit=1111,FMT='(F9.3)',ADVANCE="NO") tdot

       write(unit=1111,FMT='(F9.3)',ADVANCE="YES") ttotal
       close(unit=1111)
      endif

      call mpi_barrier(comm,infompi)

*
*     Destroy the instance of MUMPS. 
*     -----------------------------------------------------------
 444   mumps_seq%JOB = -2
      CALL DMUMPS(mumps_seq)
*
 555  call mpi_barrier(comm,infompi)
      call MPI_FINALIZE(infompi)
*
      stop


556   if(me.eq.0) print *, '########################################'
      if(me.eq.0) print *, '           ERROR Pcond solv  ',ierr
      if(me.eq.0) print *, '########################################'
      mumps_seq%JOB = -2
      CALL DMUMPS(mumps_seq)
      call mpi_barrier(comm,infompi)
      call MPI_FINALIZE(infompi)
      stop




      end
 


