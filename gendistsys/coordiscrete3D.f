       subroutine local_Discret(me,x_domains,y_domains,z_domains,
     &  size_x, size_y, size_z, a, ia, ja, nnz, p1, p2,iboundary,
     &  problem, an, ac, ad, af,right, left, top, bottom,up,down)
* Azzam modify on 13 June 2006
* input parameters
       implicit none
       integer me,x_domains,y_domains,z_domains, size_x, size_y, size_z
       integer right, left, top, bottom, up, down
       integer problem, an
       real*8  p1, p2, ac, ad, af
* output parameters
       integer ia(*),ja(*), nnz, iboundary(*)
       real*8  a(*)
* nnz                    : number of non zeros in the discretization matrix
*                          stored in coordinate format (a, ia, ja).
*
* local variables
       integer idom, idom1, size_xy, xy_domains
       integer xloc, yloc, zloc, xy_coord
       integer i, j, k, ideb, l, icr
       integer idomx, idomy, idomz
       real*8  hx,hy,hz
       real*8  s, b, c, d, e, xp, xm, yp,ym,zp,zm
       real*4  x,y,z,hx2,hy2,hz2
       integer xbg,ybg,zbg,h1x,h1y,h1z
       include 'mpif.h'

* nnz                    : number of non zeros in the discretization matrix
*                          stored in coordinate format (a, ia, ja).
* x y z                  : global space coordinate in the whole cube [0 1]
* xbg, ybg, zbg          : global starting point coordinate for each subdomain in the whole cube [0 1]
* xloc, yloc, zloc       : local coordinate according to the subdomain mesh ie according to size_x size_y size_z



*
* setup the mesh spacing, and starting point (x,y,z)
*
*      print *, me,'starting discretisation'
*      call mpi_barrier(MPI_COMM_WORLD,i)
       idom       = me+1
       xy_domains = x_domains*y_domains
       hx  = 1.0d0/dble(x_domains*(size_x-1))
       hy  = 1.0d0/dble(y_domains*(size_y-1))
       hz  = 1.0d0/dble(z_domains*(size_z-1))
       h1x = (x_domains*(size_x-1))  ! h1x=1/hx    nb of delta_x interval spacing
       h1y = (y_domains*(size_y-1))
       h1z = (z_domains*(size_z-1))
       hx2 = 0.5                     
       hy2 = 0.5                      
       hz2 = 0.5                    
       idomz = me/xy_domains
       idom1 = mod(me,xy_domains)
       idomy = idom1/x_domains
       idomx = mod(idom1,x_domains)
* compute bottom front left corner of the subdomain.
       xbg   = idomx*(size_x - 1)+1
       ybg   = idomy*(size_y - 1)+1
       zbg   = idomz*(size_z - 1)+1 
       if(me.eq.0)write(unit=6,
     &    FMT='(A,I15,I15,I15)')
     &    ' Number of nodes in xyz direction (h1x,h1y,h1z)      :',
     &    h1x,h1y,h1z
       if(me.eq.0)write(unit=6,
     &    FMT='(A,F15.3,F15.3,F15.3,I5)')
     &    ' Value of h_space in xyz direction hx,hy,hz and --p--:',
     &    hx,hy,hz,int(p1)
*
* setup the default value of iboundary (ie zero everywhere, non-boundary)
*
       do i=1,size_x*size_y*size_z
          iboundary(i) = 0
       enddo
* define the neighboring domains, and sizes
*
       size_xy   = size_x*size_y
       xy_domains = x_domains*y_domains
       right  = 0
       left   = 0
       top    = 0
       bottom = 0
       up     = 0
       down   = 0
*
       k    =  1
       if (x_domains.gt.1) then
* check if the domain has a right neighbor
         if (mod(idom,x_domains).ne.0) then
           right                  = 1
         endif
         if (mod(idom,x_domains).ne.1) then
* check if the domain has a left neighbor
           left                   = 1
         endif
       endif
*
       if (y_domains.gt.1) then
* check if the domain has a top neighbor
         idom1 = mod(idom-1,xy_domains)
         if ((idom1/x_domains).ne.(y_domains-1)) then
           top                    = 1
         endif
* check if the domain has a bottom neighbor
         if ((idom1/x_domains).ne.0) then
           bottom                 = 1
         endif
       endif
*
       if (z_domains.gt.1) then
* check if the domain has an up (above) neighbor
         if (((idom-1)/xy_domains).ne.(z_domains-1)) then
           up                     = 1
         endif
         if (((idom-1)/xy_domains).ne.0) then
* check if the domain has a down (below) neighbor
           down                   = 1
         endif
       endif
*
* setup the discretization matrix in "finite element" way
* We suppose Dirichlet boundary condition and we store the
* Dirichlet points
*
       do i=1,size_x*size_y*size_z
          iboundary(i) = 0
       enddo
*
       if (right.eq.0) then
         ideb              = size_x
         do j=1,size_z
           do i=1,size_y
             l             = ideb + (i-1)*size_x + (j-1)*size_xy
             iboundary(l)      = iboundary(l) + 1
             k = k+1
           enddo
         enddo
       endif
       if (left.eq.0) then
         ideb              = 1
         do j=1,size_z
           do i=1,size_y
             l             = ideb + (i-1)*size_x + (j-1)*size_xy
             iboundary(l)      = iboundary(l) + 1
             k             = k+1
           enddo
         enddo
       endif
*
       if (top.eq.0) then
         ideb              = (size_y-1)*size_x + 1
         do j=1,size_z
           do i=1,size_x
             l             = ideb + (j-1)*size_xy + (i-1)
             iboundary(l)      = iboundary(l) + 1
             k             = k+1
           enddo
         enddo
       endif
       if (bottom.eq.0) then
         ideb              = 1
         do j=1,size_z
           do i=1,size_x
             l             = ideb + (j-1)*size_xy + (i-1)
             iboundary(l)      = iboundary(l) + 1
             k             = k+1
           enddo
         enddo
       endif
*
       if (up.eq.0) then
         ideb              = size_xy*(size_z-1)
         do j=1,size_y
           do i=1,size_x
             l             = ideb + (j-1)*size_x + i
             iboundary(l)      = iboundary(l) + 1
             k             = k+1
           enddo
         enddo
       endif
       if (down.eq.0) then
         ideb              = 0
         do j=1,size_y
           do i=1,size_x
             l             = ideb + (j-1)*size_x + i
             iboundary(l)      = iboundary(l) + 1
             k             = k+1
           enddo
         enddo 
       endif      
*
*       s  = 6.0d0
*       b  = -(1.0d0+p1)
*       xm = -1.0d0
*       c  = -(1.0d0+p2)
*       ym = -1.0d0
*       d  = -1.0d0+p1
*       xp = -1.0d0
*       e  = -1.0d0+p2
*       yp = -1.0d0
*       zm = -1.0d0
*       zp = -1.0d0

*      Start discretization 
*      explain iboundary:
*      iboundary(i) equal the to the boundary face that lead a node
*      else for interior iboundary equal zero
      k   = 1
      icr = 1
*      if(me.eq.0)write(unit=6,
*     &        FMT='(A,I5,F14.3,F14.3,F14.3,F14.3,F14.3)')
*     &            "voici h de proc",me,hx,hy,hz,p1,p2
*      if(me.eq.1) write(unit=6,FMT='(A,I5,F14.3,F14.3,F14.3)')
*     &            "voici h de proc",me,hx,hy,hz
*
*      if(me.eq.1)write(unit=6,FMT='(A5,A5,A5,A5,A14,A14,A14)') 
*     &              'i','xloc','yloc','zloc','x','y','z'

      do i=1,size_x*size_y*size_z
         xy_coord = mod((i-1),size_xy) + 1
         xloc     = mod((xy_coord-1),size_x) + 1 ! mod((i-1),size_x) + 1 
         yloc     = ((xy_coord-1)/size_x) +1
         zloc     = ((i-1)/size_xy) + 1 
         x = xbg + xloc-1
         y = ybg + yloc-1
         z = zbg + zloc-1
         if((xloc.eq.size_x).and.(x.ne.(idomx+1)*(size_x - 1)+1))
     &      print*,'$$$$$$$$$ ERROR on discretisation 5 $$$$$$$$$'
         if((yloc.eq.size_y).and.(y.ne.(idomy+1)*(size_y - 1)+1))
     &      print*,'$$$$$$$$$ ERROR on discretisation 5 $$$$$$$$$'
         if((zloc.eq.size_z).and.(z.ne.(idomz+1)*(size_z - 1)+1))
     &      print*,'$$$$$$$$$ ERROR on discretisation 5 $$$$$$$$$'

*         if(me.eq.0)write(unit=6,FMT='(I5,I5,I5,I5,F14.3,F14.3,F14.3)') 
*     &              i,xloc,yloc,zloc,x,y,z
         if (iboundary(i).ne.0) then
*******************************************************************
********************* CASE  x = 0 BOUNDARY ************************
*******************************************************************
           if(xloc.eq.1) then
*             -----------------------------------------------------
*                                Case x=0 y=0
*             -----------------------------------------------------
              if(yloc.eq.1) then    
                 if(zloc.eq.1) then ! Corner1
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  0.25d0  
                       CASE(2) 
                               a(k)  =  0.5d0 
                       CASE(3) 
                               a(k)  =  1.0d0 
                    ENDSELECT
                 else if(zloc.eq.size_z) then ! Corner5
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  0.25d0  
                       CASE(2) 
                               a(k)  =  0.5d0 
                       CASE(3) 
                               a(k)  =  1.0d0 
                    ENDSELECT
                 else ! Edge5 
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  0.5d0  
                       CASE(2) 
                               a(k)  =  1.0d0 
                       CASE(3) 
                               stop 'ERROR ON DISCRETISATION'
                    ENDSELECT
                 endif  ! end choice on z value
*             -----------------------------------------------------
*                                Case x=0 y=size_y
*             -----------------------------------------------------
              else if(yloc.eq.size_y) then   
                 if(zloc.eq.1) then ! Corner3
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  0.25d0  
                       CASE(2) 
                               a(k)  =  0.5d0 
                       CASE(3) 
                               a(k)  =  1.0d0 
                    ENDSELECT
                 else if(zloc.eq.size_z) then ! Corner7
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  0.25d0  
                       CASE(2) 
                               a(k)  =  0.5d0 
                       CASE(3) 
                               a(k)  =  1.0d0 
                    ENDSELECT
                 else ! Edge7  
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  0.5d0  
                       CASE(2) 
                               a(k)  =  1.0d0 
                       CASE(3) 
                               stop 'ERROR ON DISCRETISATION'
                    ENDSELECT
                 endif  ! end choice on z value
*             -----------------------------------------------------
*                                Case x=0 y=else
*             -----------------------------------------------------
              else 
                 if(zloc.eq.1) then ! Edge2
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  0.5d0  
                       CASE(2) 
                               a(k)  =  1.0d0 
                       CASE(3) 
                               stop 'ERROR ON DISCRETISATION'
                    ENDSELECT
                 else if(zloc.eq.size_z) then ! Edge10
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  0.5d0  
                       CASE(2) 
                               a(k)  =  1.0d0 
                       CASE(3) 
                               stop 'ERROR ON DISCRETISATION'
                    ENDSELECT
                 else ! Face3
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  1.0d0  
                       CASE(2) 
                               stop 'ERROR ON DISCRETISATION'
                       CASE(3) 
                               stop 'ERROR ON DISCRETISATION'
                    ENDSELECT
                 endif ! end choice on z value
              endif    ! end choice on y value
*******************************************************************
****************** CASE  x = size_x BOUNDARY **********************
*******************************************************************
           else if(xloc.eq.size_x) then 
*             -----------------------------------------------------
*                                Case x=size_x y=0
*             -----------------------------------------------------
              if(yloc.eq.1) then    
                 if(zloc.eq.1) then ! Corner2
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  0.25d0  
                       CASE(2) 
                               a(k)  =  0.5d0 
                       CASE(3) 
                               a(k)  =  1.0d0 
                    ENDSELECT
                 else if(zloc.eq.size_z) then ! Corner6
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  0.25d0  
                       CASE(2) 
                               a(k)  =  0.5d0 
                       CASE(3) 
                               a(k)  =  1.0d0 
                    ENDSELECT
                 else ! Edge6
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  0.5d0  
                       CASE(2) 
                               a(k)  =  1.0d0 
                       CASE(3) 
                               stop 'ERROR ON DISCRETISATION'
                    ENDSELECT
                 endif  ! end choice on z value
*             -----------------------------------------------------
*                                Case x=size_x y=size_y
*             -----------------------------------------------------
              else if(yloc.eq.size_y) then   
                 if(zloc.eq.1) then ! Corner4
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  0.25d0  
                       CASE(2) 
                               a(k)  =  0.5d0 
                       CASE(3) 
                               a(k)  =  1.0d0 
                    ENDSELECT
                 else if(zloc.eq.size_z) then ! Corner8
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  0.25d0  
                       CASE(2) 
                               a(k)  =  0.5d0 
                       CASE(3) 
                               a(k)  =  1.0d0 
                    ENDSELECT
                 else ! Edge8
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  0.5d0  
                       CASE(2) 
                               a(k)  =  1.0d0 
                       CASE(3) 
                               stop 'ERROR ON DISCRETISATION'
                    ENDSELECT
                 endif  ! end choice on z value
*             -----------------------------------------------------
*                                Case x=size_x y=else
*             -----------------------------------------------------
              else
                 if(zloc.eq.1) then ! Edge3
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  0.5d0  
                       CASE(2) 
                               a(k)  =  1.0d0 
                       CASE(3) 
                               stop 'ERROR ON DISCRETISATION'
                    ENDSELECT
                 else if(zloc.eq.size_z) then ! Edge11
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  0.5d0  
                       CASE(2) 
                               a(k)  =  1.0d0 
                       CASE(3) 
                               stop 'ERROR ON DISCRETISATION'
                    ENDSELECT
                 else ! Face4
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  1.0d0  
                       CASE(2) 
                               stop 'ERROR ON DISCRETISATION'
                       CASE(3) 
                               stop 'ERROR ON DISCRETISATION'
                    ENDSELECT
                 endif ! end choice on z value
              endif    ! end choice on y value
*******************************************************************
******************** CASE  x = else BOUNDARY **********************
*******************************************************************
           else 
*             -----------------------------------------------------
*                               Case  x=else  y=0
*             -----------------------------------------------------
              if(yloc.eq.1) then    
                 if(zloc.eq.1) then ! Edge1
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  0.5d0  
                       CASE(2) 
                               a(k)  =  1.0d0 
                       CASE(3) 
                               stop 'ERROR ON DISCRETISATION'
                    ENDSELECT
                 else if(zloc.eq.size_z) then ! Edge9
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  0.5d0  
                       CASE(2) 
                               a(k)  =  1.0d0 
                       CASE(3) 
                               stop 'ERROR ON DISCRETISATION'
                    ENDSELECT
                 else ! Face5
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  1.0d0  
                       CASE(2) 
                               stop 'ERROR ON DISCRETISATION'
                       CASE(3) 
                               stop 'ERROR ON DISCRETISATION'
                    ENDSELECT
                 endif ! end choice on z value
*             -----------------------------------------------------
*                               Case  x=else  y=size_y
*             -----------------------------------------------------
              else if(yloc.eq.size_y) then   
                 if(zloc.eq.1) then ! Edge4
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  0.5d0  
                       CASE(2) 
                               a(k)  =  1.0d0 
                       CASE(3) 
                               stop 'ERROR ON DISCRETISATION'
                    ENDSELECT
                 else if(zloc.eq.size_z) then ! Edge12
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  0.5d0  
                       CASE(2) 
                               a(k)  =  1.0d0 
                       CASE(3) 
                               stop 'ERROR ON DISCRETISATION'
                    ENDSELECT
                 else ! Face6
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  1.0d0  
                       CASE(2) 
                               stop 'ERROR ON DISCRETISATION'
                       CASE(3) 
                               stop 'ERROR ON DISCRETISATION'
                    ENDSELECT
                 endif ! end choice on z value
*             -----------------------------------------------------
*                              Case x=else  y=else
*             -----------------------------------------------------
              else
                 if(zloc.eq.1) then ! Face1
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  1.0d0  
                       CASE(2) 
                               stop 'ERROR ON DISCRETISATION'
                       CASE(3) 
                               stop 'ERROR ON DISCRETISATION'
                    ENDSELECT
                 else if(zloc.eq.size_z) then ! Face2
                    SELECTCASE(iboundary(i))
                       CASE(1) 
                               a(k)  =  1.0d0  
                       CASE(2) 
                               stop 'ERROR ON DISCRETISATION'
                       CASE(3) 
                               stop 'ERROR ON DISCRETISATION'
                    ENDSELECT
                 else ! Interior IMPOSSIBLE D'ENTRER DANS CETTE CONDITION
                    stop'ERROR ON DISCRETISATION'
                 endif ! end choice on z value
              endif ! end choice on y value
           endif ! end choice on x value
           ia(k) = icr
           ja(k) = icr
           k     = k+1
           icr   = icr+1
         else
*******************************************************************
************************* CASE  x = 0 *****************************
*******************************************************************
           if(xloc.eq.1) then
*             -----------------------------------------------------
*                                Case x=0 y=0
*             -----------------------------------------------------
              if(yloc.eq.1) then    
                 if(zloc.eq.1) then ! Corner1
                    call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                 s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                      hx,hy,hz)
                    a(k)  = s/8.0d0         ! diag coef
                    ia(k) = icr
                    ja(k) = icr
                    k     = k+1
                    a(k)  = zp/4.0d0        ! z-up point zp
                    ia(k) = icr
                    ja(k) = icr + size_xy
                    k     = k+1
                    a(k)  = yp/4.0d0        ! north point yp
                    ia(k) = icr
                    ja(k) = icr + size_x
                    k     = k+1
                    a(k)  = xp/4.0d0        ! east point xp
                    ia(k) = icr
                    ja(k) = icr + 1
                    k     = k+1
                    icr   = icr+1 !finish this node
                 else
                    if(zloc.eq.size_z) then ! Corner5
                       call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                    s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                         hx,hy,hz)

                       a(k)  = s/8.0d0         ! diag coef
                       ia(k) = icr
                       ja(k) = icr
                       k     = k+1
                       a(k)  = zm/4.0d0        ! z-down point zm
                       ia(k) = icr
                       ja(k) = icr - size_xy
                       k     = k+1
                       a(k)  = yp/4.0d0        ! north point yp
                       ia(k) = icr
                       ja(k) = icr + size_x
                       k     = k+1
                       a(k)  = xp/4.0d0        ! east point xp
                       ia(k) = icr
                       ja(k) = icr + 1
                       k     = k+1
                       icr   = icr+1  !finish this node
                    else ! Edge5 ! Je teste le iboundary seulement pour direction z pour savoir si connecter a un noeud boundary et de meme pour les autres
                       call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                    s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                         hx,hy,hz)

                       a(k)  = s/4.0d0         ! diag coef
                       ia(k) = icr
                       ja(k) = icr
                       k     = k+1
                       a(k)  = zp/4.0d0        ! z-up point zp
                       ia(k) = icr
                       ja(k) = icr + size_xy
                       if(iboundary(ja(k)).eq.0) k     = k+1      
                       a(k)  = zm/4.0d0        ! z-down point zm
                       ia(k) = icr
                       ja(k) = icr - size_xy
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = yp/2.0d0        ! north point yp
                       ia(k) = icr
                       ja(k) = icr + size_x
                       k     = k+1
                       a(k)  = xp/2.0d0        ! east point xp
                       ia(k) = icr
                       ja(k) = icr + 1
                       k     = k+1
                       icr   = icr+1  !finish this node
                    endif
                 endif  ! end choice on z value
*             -----------------------------------------------------
*                                Case x=0 y=size_y
*             -----------------------------------------------------
              else
              if(yloc.eq.size_y) then   
                 if(zloc.eq.1) then ! Corner3
                    call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                 s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                      hx,hy,hz)

                    a(k)  = s/8.0d0         ! diag coef
                    ia(k) = icr
                    ja(k) = icr
                    k     = k+1
                    a(k)  = zp/4.0d0        ! z-up point zp
                    ia(k) = icr
                    ja(k) = icr + size_xy
                    k     = k+1
                    a(k)  = ym/4.0d0        ! south point  ym
                    ia(k) = icr
                    ja(k) = icr - size_x
                    k     = k+1
                    a(k)  = xp/4.0d0        ! east point xp
                    ia(k) = icr
                    ja(k) = icr + 1
                    k     = k+1
                    icr   = icr+1 !finish this node
                 else
                    if(zloc.eq.size_z) then ! Corner7
                       call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                    s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                         hx,hy,hz)

                       a(k)  = s/8.0d0         ! diag coef
                       ia(k) = icr
                       ja(k) = icr
                       k     = k+1
                       a(k)  = zm/4.0d0        ! z-down point zm
                       ia(k) = icr
                       ja(k) = icr - size_xy
                       k     = k+1
                       a(k)  = ym/4.0d0        ! south point  ym
                       ia(k) = icr
                       ja(k) = icr - size_x
                       k     = k+1
                       a(k)  = xp/4.0d0        ! east point xp
                       ia(k) = icr
                       ja(k) = icr + 1
                       k     = k+1
                       icr   = icr+1  !finish this node
                    else ! Edge7 ! Je teste le iboundary seulement pour direction z pour savoir si connecter a un noeud boundary et de meme pour les autres 
                       call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                    s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                         hx,hy,hz)

                       a(k)  = s/4.0d0         ! diag coef
                       ia(k) = icr
                       ja(k) = icr
                       k     = k+1
                       a(k)  = zp/4.0d0        ! z-up point zp
                       ia(k) = icr
                       ja(k) = icr + size_xy
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = zm/4.0d0        ! z-down point zm
                       ia(k) = icr
                       ja(k) = icr - size_xy
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = ym/2.0d0        ! south point  ym
                       ia(k) = icr
                       ja(k) = icr - size_x
                       k     = k+1
                       a(k)  = xp/2.0d0        ! east point xp
                       ia(k) = icr
                       ja(k) = icr + 1
                       k     = k+1
                       icr   = icr+1 !finish this node
                    endif
                 endif  ! end choice on z value
*             -----------------------------------------------------
*                                Case x=0 y=else
*             -----------------------------------------------------
              else
                 if(zloc.eq.1) then ! Edge2
                    call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                 s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                      hx,hy,hz)

                    a(k)  = s/4.0d0         ! diag coef
                    ia(k) = icr
                    ja(k) = icr
                    k     = k+1
                    a(k)  = zp/2.0d0        ! z-up point zp
                    ia(k) = icr
                    ja(k) = icr + size_xy
                    k     = k+1
                    a(k)  = yp/4.0d0        ! north point yp
                    ia(k) = icr
                    ja(k) = icr + size_x
                    if(iboundary(ja(k)).eq.0) k     = k+1
                    a(k)  = ym/4.0d0        ! south point  ym
                    ia(k) = icr
                    ja(k) = icr - size_x
                    if(iboundary(ja(k)).eq.0) k     = k+1
                    a(k)  = xp/2.0d0        ! east point xp
                    ia(k) = icr
                    ja(k) = icr + 1
                    k     = k+1
                    icr   = icr+1 !finish this node
                 else
                    if(zloc.eq.size_z) then ! Edge10
                       call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                    s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                         hx,hy,hz)
                       a(k)  = s/4.0d0         ! diag coef
                       ia(k) = icr
                       ja(k) = icr
                       k     = k+1
                       a(k)  = zm/2.0d0        ! z-down point zm
                       ia(k) = icr
                       ja(k) = icr - size_xy
                       k     = k+1
                       a(k)  = yp/4.0d0        ! north point yp
                       ia(k) = icr
                       ja(k) = icr + size_x
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = ym/4.0d0        ! south point  ym
                       ia(k) = icr
                       ja(k) = icr - size_x
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = xp/2.0d0        ! east point xp
                       ia(k) = icr
                       ja(k) = icr + 1
                       k     = k+1
                       icr   = icr+1  !finish this node
                    else ! Face3
                       call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                    s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                         hx,hy,hz)

                       a(k)  = s/2.0d0         ! diag coef
                       ia(k) = icr
                       ja(k) = icr
                       k     = k+1
                       a(k)  = zp/2.0d0        ! z-up point zp
                       ia(k) = icr
                       ja(k) = icr + size_xy
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = zm/2.0d0        ! z-down point zm
                       ia(k) = icr
                       ja(k) = icr - size_xy
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = yp/2.0d0        ! north point yp
                       ia(k) = icr
                       ja(k) = icr + size_x
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = ym/2.0d0        ! south point  ym
                       ia(k) = icr
                       ja(k) = icr - size_x
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = xp             ! east point xp  ! interior node
                       ia(k) = icr
                       ja(k) = icr + 1
                       k     = k+1
                       icr   = icr+1  !finish this node
                    endif
                 endif ! end choice on z value
              endif
              endif    ! end choice on y value
*******************************************************************
*********************** CASE  x = size_x **************************
*******************************************************************
           else
           if(xloc.eq.size_x) then 
*             -----------------------------------------------------
*                                Case x=size_x y=0
*             -----------------------------------------------------
              if(yloc.eq.1) then    
                 if(zloc.eq.1) then ! Corner2
                    call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                 s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                      hx,hy,hz)

                    a(k)  = s/8.0d0         ! diag coef
                    ia(k) = icr
                    ja(k) = icr
                    k     = k+1
                    a(k)  = zp/4.0d0        ! z-up point zp
                    ia(k) = icr
                    ja(k) = icr + size_xy
                    k     = k+1
                    a(k)  = yp/4.0d0        ! north point yp
                    ia(k) = icr
                    ja(k) = icr + size_x
                    k     = k+1
                    a(k)  = xm/4.0d0        ! west point xm
                    ia(k) = icr
                    ja(k) = icr - 1
                    k     = k+1
                    icr   = icr+1 !finish this node
                 else
                    if(zloc.eq.size_z) then ! Corner6
                       call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                    s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                         hx,hy,hz)

                       a(k)  = s/8.0d0         ! diag coef
                       ia(k) = icr
                       ja(k) = icr
                       k     = k+1
                       a(k)  = zm/4.0d0        ! z-down point zm
                       ia(k) = icr
                       ja(k) = icr - size_xy
                       k     = k+1
                       a(k)  = yp/4.0d0        ! north point yp
                       ia(k) = icr
                       ja(k) = icr + size_x
                       k     = k+1
                       a(k)  = xm/4.0d0        ! west point xm
                       ia(k) = icr
                       ja(k) = icr - 1
                       k     = k+1
                       icr   = icr+1 !finish this node
                    else ! Edge6 
                       call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                    s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                         hx,hy,hz)

                       a(k)  = s/4.0d0         ! diag coef
                       ia(k) = icr
                       ja(k) = icr
                       k     = k+1
                       a(k)  = zp/4.0d0        ! z-up point zp
                       ia(k) = icr
                       ja(k) = icr + size_xy
                       if(iboundary(ja(k)).eq.0) k     = k+1      
                       a(k)  = zm/4.0d0        ! z-down point zm
                       ia(k) = icr
                       ja(k) = icr - size_xy
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = yp/2.0d0        ! north point yp
                       ia(k) = icr
                       ja(k) = icr + size_x
                       k     = k+1
                       a(k)  = xm/2.0d0        ! west point xm
                       ia(k) = icr
                       ja(k) = icr - 1
                       k     = k+1
                       icr   = icr+1 !finish this node
                    endif
                 endif  ! end choice on z value
*             -----------------------------------------------------
*                                Case x=size_x y=size_y
*             -----------------------------------------------------
              else
              if(yloc.eq.size_y) then   
                 if(zloc.eq.1) then ! Corner4
                    call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                  s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                       hx,hy,hz)

                    a(k)  = s/8.0d0         ! diag coef
                    ia(k) = icr
                    ja(k) = icr
                    k     = k+1
                    a(k)  = zp/4.0d0        ! z-up point zp
                    ia(k) = icr
                    ja(k) = icr + size_xy
                    k     = k+1
                    a(k)  = ym/4.0d0        ! south point  ym
                    ia(k) = icr
                    ja(k) = icr - size_x
                    k     = k+1
                    a(k)  = xm/4.0d0        ! west point xm
                    ia(k) = icr
                    ja(k) = icr - 1
                    k     = k+1
                    icr   = icr+1 !finish this node
                 else
                    if(zloc.eq.size_z) then ! Corner8
                       call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                    s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                         hx,hy,hz)

                       a(k)  = s/8.0d0         ! diag coef
                       ia(k) = icr
                       ja(k) = icr
                       k     = k+1
                       a(k)  = zm/4.0d0        ! z-down point zm
                       ia(k) = icr
                       ja(k) = icr - size_xy
                       k     = k+1
                       a(k)  = ym/4.0d0        ! south point  ym
                       ia(k) = icr
                       ja(k) = icr - size_x
                       k     = k+1
                       a(k)  = xm/4.0d0        ! west point xm
                       ia(k) = icr
                       ja(k) = icr - 1
                       k     = k+1
                       icr   = icr+1 !finish this node
                    else ! Edge8 
                       call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                    s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                         hx,hy,hz)

                       a(k)  = s/4.0d0         ! diag coef
                       ia(k) = icr
                       ja(k) = icr
                       k     = k+1
                       a(k)  = zp/4.0d0        ! z-up point zp
                       ia(k) = icr
                       ja(k) = icr + size_xy
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = zm/4.0d0        ! z-down point zm
                       ia(k) = icr
                       ja(k) = icr - size_xy
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = ym/2.0d0        ! south point  ym
                       ia(k) = icr
                       ja(k) = icr - size_x
                       k     = k+1
                       a(k)  = xm/2.0d0        ! west point xm
                       ia(k) = icr
                       ja(k) = icr - 1
                       k     = k+1
                       icr   = icr+1 !finish this node
                    endif
                 endif  ! end choice on z value
*             -----------------------------------------------------
*                                Case x=size_x y=else
*             -----------------------------------------------------
              else
                 if(zloc.eq.1) then ! Edge3
                    call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                 s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                      hx,hy,hz)

                    a(k)  = s/4.0d0         ! diag coef
                    ia(k) = icr
                    ja(k) = icr
                    k     = k+1
                    a(k)  = zp/2.0d0        ! z-up point zp
                    ia(k) = icr
                    ja(k) = icr + size_xy
                    k     = k+1
                    a(k)  = yp/4.0d0        ! north point yp
                    ia(k) = icr
                    ja(k) = icr + size_x
                    if(iboundary(ja(k)).eq.0) k     = k+1
                    a(k)  = ym/4.0d0        ! south point  ym
                    ia(k) = icr
                    ja(k) = icr - size_x
                    if(iboundary(ja(k)).eq.0) k     = k+1
                    a(k)  = xm/2.0d0        ! west point xm
                    ia(k) = icr
                    ja(k) = icr - 1
                    k     = k+1
                    icr   = icr+1 !finish this node
                 else
                    if(zloc.eq.size_z) then ! Edge11
                       call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                    s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                         hx,hy,hz)

                       a(k)  = s/4.0d0         ! diag coef
                       ia(k) = icr
                       ja(k) = icr
                       k     = k+1
                       a(k)  = zm/2.0d0        ! z-down point zm
                       ia(k) = icr
                       ja(k) = icr - size_xy
                       k     = k+1
                       a(k)  = yp/4.0d0        ! north point yp
                       ia(k) = icr
                       ja(k) = icr + size_x
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = ym/4.0d0        ! south point  ym
                       ia(k) = icr
                       ja(k) = icr - size_x
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = xm/2.0d0        ! west point xm
                       ia(k) = icr
                       ja(k) = icr - 1
                       k     = k+1
                       icr   = icr+1 !finish this node
                    else ! Face4
                       call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                    s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                         hx,hy,hz)

                       a(k)  = s/2.0d0         ! diag coef
                       ia(k) = icr
                       ja(k) = icr
                       k     = k+1
                       a(k)  = zp/2.0d0        ! z-up point zp
                       ia(k) = icr
                       ja(k) = icr + size_xy
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = zm/2.0d0        ! z-down point zm
                       ia(k) = icr
                       ja(k) = icr - size_xy
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = yp/2.0d0        ! north point yp
                       ia(k) = icr
                       ja(k) = icr + size_x
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = ym/2.0d0        ! south point  ym
                       ia(k) = icr
                       ja(k) = icr - size_x
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = xm             ! west point xm  ! interior node
                       ia(k) = icr
                       ja(k) = icr - 1
                       k     = k+1
                       icr   = icr+1 !finish this node
                    endif
                 endif ! end choice on z value
              endif
              endif    ! end choice on y value
*******************************************************************
****************** CASE  x = else (interior x) ********************
*******************************************************************
           else 
*             -----------------------------------------------------
*                               Case  x=else  y=0
*             -----------------------------------------------------
              if(yloc.eq.1) then    
                 if(zloc.eq.1) then ! Edge1
                    call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                 s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                      hx,hy,hz)

                    a(k)  = s/4.0d0         ! diag coef
                    ia(k) = icr
                    ja(k) = icr
                    k     = k+1
                    a(k)  = zp/2.0d0        ! z-up point zp
                    ia(k) = icr
                    ja(k) = icr + size_xy
                    k     = k+1
                    a(k)  = yp/2.0d0        ! north point yp
                    ia(k) = icr
                    ja(k) = icr + size_x
                    k     = k+1
                    a(k)  = xp/4.0d0        ! east point xp
                    ia(k) = icr
                    ja(k) = icr + 1
                    if(iboundary(ja(k)).eq.0) k     = k+1
                    a(k)  = xm/4.0d0        ! west point xm  
                    ia(k) = icr
                    ja(k) = icr - 1
                    if(iboundary(ja(k)).eq.0) k     = k+1
                    icr   = icr+1 !finish this node
                 else
                    if(zloc.eq.size_z) then ! Edge9
                       call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                    s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                         hx,hy,hz)

                       a(k)  = s/4.0d0         ! diag coef
                       ia(k) = icr
                       ja(k) = icr
                       k     = k+1
                       a(k)  = zm/2.0d0        ! z-down point zm
                       ia(k) = icr
                       ja(k) = icr - size_xy
                       k     = k+1
                       a(k)  = yp/2.0d0        ! north point yp
                       ia(k) = icr
                       ja(k) = icr + size_x
                       k     = k+1
                       a(k)  = xp/4.0d0        ! east point xp
                       ia(k) = icr
                       ja(k) = icr + 1
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = xm/4.0d0        ! west point xm  
                       ia(k) = icr
                       ja(k) = icr - 1
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       icr   = icr+1 !finish this node
                    else ! Face5
                       call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                    s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                         hx,hy,hz)

                       a(k)  = s/2.0d0         ! diag coef
                       ia(k) = icr
                       ja(k) = icr
                       k     = k+1
                       a(k)  = zp/2.0d0        ! z-up point zp
                       ia(k) = icr
                       ja(k) = icr + size_xy
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = zm/2.0d0        ! z-down point zm
                       ia(k) = icr
                       ja(k) = icr - size_xy
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = yp              ! north point yp ! interior node
                       ia(k) = icr
                       ja(k) = icr + size_x
                       k     = k+1
                       a(k)  = xp/2.0d0        ! east point xp  
                       ia(k) = icr
                       ja(k) = icr + 1
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = xm/2.0d0        ! west point xm  
                       ia(k) = icr
                       ja(k) = icr - 1
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       icr   = icr+1  !finish this node
                    endif
                 endif ! end choice on z value
*             -----------------------------------------------------
*                               Case  x=else  y=size_y
*             -----------------------------------------------------
              else
              if(yloc.eq.size_y) then   
                 if(zloc.eq.1) then ! Edge4
                    call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                 s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                      hx,hy,hz)

                    a(k)  = s/4.0d0         ! diag coef
                    ia(k) = icr
                    ja(k) = icr
                    k     = k+1
                    a(k)  = zp/2.0d0        ! z-up point zp
                    ia(k) = icr
                    ja(k) = icr + size_xy
                    k     = k+1
                    a(k)  = ym/2.0d0        ! south point ym
                    ia(k) = icr
                    ja(k) = icr - size_x
                    k     = k+1
                    a(k)  = xp/4.0d0        ! east point xp
                    ia(k) = icr
                    ja(k) = icr + 1
                    if(iboundary(ja(k)).eq.0) k     = k+1
                    a(k)  = xm/4.0d0        ! west point xm  
                    ia(k) = icr
                    ja(k) = icr - 1
                    if(iboundary(ja(k)).eq.0) k     = k+1
                    icr   = icr+1 !finish this node
                 else
                    if(zloc.eq.size_z) then ! Edge12
                       call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                    s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                         hx,hy,hz)

                       a(k)  = s/4.0d0         ! diag coef
                       ia(k) = icr
                       ja(k) = icr
                       k     = k+1
                       a(k)  = zm/2.0d0        ! z-down point zm
                       ia(k) = icr
                       ja(k) = icr - size_xy
                       k     = k+1
                       a(k)  = ym/2.0d0        ! south point ym
                       ia(k) = icr
                       ja(k) = icr - size_x
                       k     = k+1
                       a(k)  = xp/4.0d0        ! east point xp
                       ia(k) = icr
                       ja(k) = icr + 1
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = xm/4.0d0        ! west point xm  
                       ia(k) = icr
                       ja(k) = icr - 1
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       icr   = icr+1 !finish this node
                    else ! Face6
                       call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                    s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                         hx,hy,hz)

                       a(k)  = s/2.0d0         ! diag coef
                       ia(k) = icr
                       ja(k) = icr
                       k     = k+1
                       a(k)  = zp/2.0d0        ! z-up point zp
                       ia(k) = icr
                       ja(k) = icr + size_xy
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = zm/2.0d0        ! z-down point zm
                       ia(k) = icr
                       ja(k) = icr - size_xy
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = ym              ! south point ym  ! interior node
                       ia(k) = icr
                       ja(k) = icr - size_x
                       k     = k+1
                       a(k)  = xp/2.0d0        ! east point xp  
                       ia(k) = icr
                       ja(k) = icr + 1
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = xm/2.0d0        ! west point xm  
                       ia(k) = icr
                       ja(k) = icr - 1
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       icr   = icr+1  !finish this node
                    endif
                 endif ! end choice on z value
*             -----------------------------------------------------
*                              Case x=else  y=else
*             -----------------------------------------------------
              else
                 if(zloc.eq.1) then ! Face1
                    call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                 s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                      hx,hy,hz)

                    a(k)  = s/2.0d0         ! diag coef
                    ia(k) = icr
                    ja(k) = icr
                    k     = k+1
                    a(k)  = zp              ! z-up point zp  ! interior node
                    ia(k) = icr
                    ja(k) = icr + size_xy
                    k     = k+1
                    a(k)  = yp/2.0d0        ! north point yp
                    ia(k) = icr
                    ja(k) = icr + size_x
                    if(iboundary(ja(k)).eq.0) k     = k+1
                    a(k)  = ym/2.0d0        ! south point  ym
                    ia(k) = icr
                    ja(k) = icr - size_x
                    if(iboundary(ja(k)).eq.0) k     = k+1
                    a(k)  = xp/2.0d0        ! east point xp  
                    ia(k) = icr
                    ja(k) = icr + 1
                    if(iboundary(ja(k)).eq.0) k     = k+1
                    a(k)  = xm/2.0d0        ! west point xm  
                    ia(k) = icr
                    ja(k) = icr - 1
                    if(iboundary(ja(k)).eq.0) k     = k+1
                    icr   = icr+1  !finish this node
                 else
                    if(zloc.eq.size_z) then ! Face2
                       call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                    s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                         hx,hy,hz)
                       a(k)  = s/2.0d0         ! diag coef
                       ia(k) = icr
                       ja(k) = icr
                       k     = k+1
                       a(k)  = zm              ! z-down point zm ! interior node
                       ia(k) = icr
                       ja(k) = icr - size_xy
                       k     = k+1
                       a(k)  = yp/2.0d0        ! north point yp
                       ia(k) = icr
                       ja(k) = icr + size_x
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = ym/2.0d0        ! south point  ym
                       ia(k) = icr
                       ja(k) = icr - size_x
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = xp/2.0d0        ! east point xp  
                       ia(k) = icr
                       ja(k) = icr + 1
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = xm/2.0d0        ! west point xm  
                       ia(k) = icr
                       ja(k) = icr - 1
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       icr   = icr+1  !finish this node
                    else ! Interior
                       call coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &                    s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                                         hx,hy,hz)

                       a(k)  = s               ! diag coef
                       ia(k) = icr
                       ja(k) = icr
                       k     = k+1
                       a(k)  = zp              ! z-up point zp   ! interior node
                       ia(k) = icr
                       ja(k) = icr + size_xy
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = zm              ! z-down point zm ! interior node
                       ia(k) = icr
                       ja(k) = icr - size_xy
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = yp              ! north point yp  ! interior node
                       ia(k) = icr
                       ja(k) = icr + size_x
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = ym              ! south point ym  ! interior node
                       ia(k) = icr
                       ja(k) = icr - size_x
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = xp              ! east point xp   ! interior node
                       ia(k) = icr
                       ja(k) = icr + 1
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       a(k)  = xm             ! west point xm    ! interior node
                       ia(k) = icr
                       ja(k) = icr - 1
                       if(iboundary(ja(k)).eq.0) k     = k+1
                       icr   = icr+1  !finish this node
                    endif
                 endif ! end choice on z value
              endif
              endif ! end choice on y value
           endif
           endif ! end choice on x value
         endif ! end if on boundary (0 or not)
      enddo ! END LOOP on iboundary(n)




*
       nnz = k - 1
*
       return
       end
*

