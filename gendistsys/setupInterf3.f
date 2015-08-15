C>
C! @param [out] sizeinterf_lo
C!   number of nodes on the interface which are handled by this domain.
C! @param [out] list_Intrfc_lo
C!   list of nodes on the interface which are handled by this domain.
C!   It must be allocated by the caller to a minimum size "sizeIntrf"
C! @param [in,out] iwork
C!   Integer workspace of size size_x * size_y * size_z 
C!
       subroutine setupInterf(me,x_domains,y_domains,z_domains,size_x,
     &       size_y,size_z,nbvi,indexVi,ptr_Index_Intrfc, Index_Intrfc,
     &       list_Intrfc,SizeInterf_lo, list_Intrfc_lo,
     &       perm, invPerm, weight, iwork)
*
* input parameters
       integer me,x_domains,y_domains,z_domains,size_x, size_y,size_z
* output parameters
       integer nbvi, indexVi(*), ptr_Index_Intrfc(*), Index_Intrfc(*)
       integer list_Intrfc(*), perm(*), invPerm(*), iwork(*)
       integer SizeInterf_lo, list_Intrfc_lo(*)
       real*8  weight(*)
* nbvi                   : number of subdomains neighbor of the subdomain 
*                          no (me+1) as me is the MPI process number starting
*                          at 0.
* indexVi(nbvi)          : no of the neighboring subdomain
* ptr_Index_Intrfc(nbvi) : pointer to the index of the first point on the
*                          interface with the sudomain. The indices of the
*                          interface points are stored in the array Index_Intrfc
*                          The indices are given in the original ordering.
*
* Index_Intrfc(*)        : list of the interface points, pointed by the
*                          array ptr_Index_Intrfc.
*                          The indices are given in the original ordering.
*
* list_Intrfc(sizeIntrf) : indices of the interface points.
*                          The indices are given in the original ordering.
*
* perm(*)                : permutation allowing to go from the "schur ordering"
*                          to the original ordering
*
* invPerm(*)             : permutation allowing to go from the original ordering
*                          to the Schur ordering where the interface points
*                          are stored first
*
* weight(sizeIntrf)      : weight given to each interface point to be used
*                          in the dot product computation. The weight is
*                          equal to the inverse of number of subdomains to
*                          which the point belongs.
*
* local variable
       integer sizeIntrf, idom, idom1, k2
       integer xy_domains, size_xy,ios
       integer i, j ,k , ideb, l, right, left, top, bottom, up, down,icr
       integer ,dimension(:), allocatable :: iwork2
       Integer :: tmp
*
       ALLOCATE(iwork2(size_x*size_y*size_z))
*
* setup the default value of weight (ie one everywhere)
*
       idom = me+1
*
       do i=1,size_x*size_y*size_z
          iwork(i) = 0
          iwork2(i) = -1
       enddo
* define the neighboring domains, and sizes
*
       size_xy    =  size_x*size_y
       xy_domains = x_domains*y_domains
       right  = 0
       left   = 0
       top    = 0
       bottom = 0
       up     = 0
       down   = 0
*
       nbvi = 0
       k    = 1
*
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
* Store the neighbors in an increasing rank order
* (will be used to generate the global ordering of the
* complete distributed matrix for MUMPS)
*
* right, left (x-dir); top, bottom (y-dir); up, down (z-dir)
       k2   =  1
*
*****************************************
* Do the upper (z) 9 neighbour domains
*****************************************
       if ((up.eq.1).and.(right.eq.1).and.(top.eq.1)) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me + xy_domains + x_domains + 1
         ptr_Index_Intrfc(nbvi) = k
         l                      = size_z*size_xy
         Index_Intrfc(k)        = l
         k                      = k + 1
         iwork2(l)              = me + xy_domains + x_domains + 1
       endif
*
       if ((up.eq.1).and.(top.eq.1)) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me + xy_domains + x_domains
         ptr_Index_Intrfc(nbvi) = k
         ideb                   = size_z*size_xy - size_x + 1
         do i = 1,size_x
           l                    = ideb + (i-1)
           Index_Intrfc(k)      = l
           k                    = k + 1
           iwork2(l)            = me + xy_domains + x_domains
         enddo
       endif
*  
       if ((up.eq.1).and.(left.eq.1).and.(top.eq.1)) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me + xy_domains + x_domains - 1
         ptr_Index_Intrfc(nbvi) = k
         l                      = size_z*size_xy - size_x + 1
         Index_Intrfc(k)        = l
         k                      = k + 1
         iwork2(l)              = me + xy_domains + x_domains - 1
       endif  
*
       if ((up.eq.1).and.(right.eq.1)) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me + xy_domains + 1
         ptr_Index_Intrfc(nbvi) = k
         ideb                   = (size_z-1)*size_xy + size_x
         do i = 1,size_y
           l                    = ideb + (i-1)*size_x
           Index_Intrfc(k)      = l
           k                    = k + 1
           iwork2(l)            = me + xy_domains + 1
         enddo
       endif  
*
       if (up.eq.1) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me + xy_domains
         ptr_Index_Intrfc(nbvi) = k
         ideb                   = (size_z-1)*size_xy + 1
         do j = 1,size_y
           do i = 1,size_x
             l                  = ideb + (i-1) + (j-1)*size_x
             Index_Intrfc(k)    = l
             k                  = k + 1
             iwork2(l)          = me + xy_domains
             iwork(l)           = iwork(l) + 1
           enddo
         enddo
       endif
*
       if ((up.eq.1).and.(left.eq.1)) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me + xy_domains - 1
         ptr_Index_Intrfc(nbvi) = k
         ideb                   = (size_z-1)*size_xy + 1
         do i = 1,size_y
           l                    = ideb + (i-1)*size_x
           Index_Intrfc(k)      = l
           k                    = k + 1
           iwork2(l)            = me + xy_domains - 1
         enddo
       endif
*
       if ((up.eq.1).and.(right.eq.1).and.(bottom.eq.1)) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me + xy_domains - x_domains + 1
         ptr_Index_Intrfc(nbvi) = k
         l                      = (size_z-1)*size_xy + size_x
         Index_Intrfc(k)        = l
         k                      = k + 1
         iwork2(l)              = me + xy_domains - x_domains + 1
       endif
*
       if ((up.eq.1).and.(bottom.eq.1)) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me + xy_domains - x_domains
         ptr_Index_Intrfc(nbvi) = k
         ideb                   = (size_z-1)*size_xy + 1
         do i = 1,size_x
           l                    = ideb + (i-1)
           Index_Intrfc(k)      = l
           k                    = k + 1
           iwork2(l)            = me + xy_domains - x_domains
         enddo
       endif 
*
       if ((up.eq.1).and.(left.eq.1).and.(bottom.eq.1)) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me + xy_domains - x_domains - 1
         ptr_Index_Intrfc(nbvi) = k
         l                      = (size_z-1)*size_xy + 1
         Index_Intrfc(k)        = l
         k                      = k + 1
         iwork2(l)              = me + xy_domains - x_domains - 1
       endif           
****************************************************************
* Do the 8 central neighbours shared by the domain- same z-level
****************************************************************
       if ((top.eq.1).and.(right.eq.1)) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me + x_domains + 1
         ptr_Index_Intrfc(nbvi) = k
         ideb                   = size_x*size_y
         do j = 1,size_z
           l                    = ideb + (j-1)*size_xy
           Index_Intrfc(k)      = l
           k                    = k + 1
           iwork2(l)            = me + x_domains + 1
         enddo
       endif
*
       if (top.eq.1) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me + x_domains
         ptr_Index_Intrfc(nbvi) = k
         ideb                   = (size_y-1)*size_x + 1
         do j = 1,size_z
           do i = 1,size_x
             l                  = ideb + (i-1) + (j-1)*size_xy
             Index_Intrfc(k)    = l
             k                  = k + 1
             iwork2(l)          = me + x_domains
             iwork(l)           = iwork(l) + 1
           enddo
         enddo
       endif
*
       if ((top.eq.1).and.(left.eq.1)) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me + x_domains - 1
         ptr_Index_Intrfc(nbvi) = k
         ideb                   = (size_y-1)*size_x + 1
         do j = 1,size_z
           l                    = ideb + (j-1)*size_xy
           Index_Intrfc(k)      = l
           k                    = k + 1
           iwork2(l)            = me + x_domains - 1
         enddo
       endif
*
       if (right.eq.1) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me + 1
         ptr_Index_Intrfc(nbvi) = k
         ideb                   = size_x
         do j = 1,size_z
           do i = 1,size_y
             l                  = ideb + (i-1)*size_x + (j-1)*size_xy
             Index_Intrfc(k)    = l
             k                  = k + 1
             iwork2(l)          = me + 1
             iwork(l)           = iwork(l) + 1
           enddo
         enddo
       endif
*
       if (left.eq.1) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me - 1
         ptr_Index_Intrfc(nbvi) = k
         ideb                   = 1
         do j = 1,size_z
           do i = 1,size_y
             l                  = ideb + (i-1)*size_x + (j-1)*size_xy
             Index_Intrfc(k)    = l
             k                  = k + 1
             iwork2(l)          = me - 1
             iwork(l)           = iwork(l) + 1
           enddo
         enddo
       endif
*
       if ((bottom.eq.1).and.(right.eq.1)) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me - x_domains + 1
         ptr_Index_Intrfc(nbvi) = k
         ideb                   = size_x
         do j = 1,size_z
           l                    = ideb + (j-1)*size_xy
           Index_Intrfc(k)      = l
           k                    = k + 1
           iwork2(l)            = me - x_domains + 1
         enddo
       endif
*
       if (bottom.eq.1) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me - x_domains
         ptr_Index_Intrfc(nbvi) = k
         ideb                   = 1
         do j = 1,size_z
           do i = 1,size_x
             l                  = ideb + (i-1) + (j-1)*size_xy
             Index_Intrfc(k)    = l
             k                  = k + 1
             iwork2(l)          = me - x_domains
             iwork(l)           = iwork(l) + 1
           enddo
         enddo
       endif
*
       if ((bottom.eq.1).and.(left.eq.1)) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me - x_domains - 1
         ptr_Index_Intrfc(nbvi) = k
         ideb                   = 1
         do j = 1,size_z
           l                    = ideb + (j-1)*size_xy
           Index_Intrfc(k)      = l
           k                    = k + 1
           iwork2(l)            = me - x_domains - 1
         enddo
       endif
*****************************************
* Do the lower (z) 9 neighbour domains
*****************************************
       if ((down.eq.1).and.(right.eq.1).and.(top.eq.1)) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me - xy_domains + x_domains + 1
         ptr_Index_Intrfc(nbvi) = k
         l                      = size_xy
         Index_Intrfc(k)        = l
         k                      = k + 1
         iwork2(l)              = me - xy_domains + x_domains + 1
       endif
*
       if ((down.eq.1).and.(top.eq.1)) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me - xy_domains + x_domains
         ptr_Index_Intrfc(nbvi) = k
         ideb                   = (size_y - 1)*size_x + 1
         do i = 1,size_x
           l                    = ideb + (i-1)
           Index_Intrfc(k)      = l
           k                    = k + 1
           iwork2(l)            = me - xy_domains + x_domains
         enddo
       endif
*
       if ((down.eq.1).and.(left.eq.1).and.(top.eq.1)) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me - xy_domains + x_domains - 1
         ptr_Index_Intrfc(nbvi) = k
         l                      = size_xy - size_x + 1
         Index_Intrfc(k)        = l
         k                      = k + 1
         iwork2(l)              = me - xy_domains + x_domains - 1
       endif
*      
       if ((down.eq.1).and.(right.eq.1)) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me - xy_domains + 1
         ptr_Index_Intrfc(nbvi) = k
         ideb                   = size_x
         do i = 1,size_y
           l                    = ideb + (i-1)*size_x
           Index_Intrfc(k)      = l
           k                    = k + 1
           iwork2(l)            = me - xy_domains + 1
         enddo
       endif  
*
       if (down.eq.1) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me - xy_domains
         ptr_Index_Intrfc(nbvi) = k
         ideb                   = 1
         do j = 1,size_y
           do i = 1,size_x
             l                  = ideb + (i-1) + (j-1)*size_x
             Index_Intrfc(k)    = l
             k                  = k + 1
             iwork2(l)          = me - xy_domains
             iwork(l)           = iwork(l) + 1
           enddo
         enddo
       endif
*
       if ((down.eq.1).and.(left.eq.1)) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me - xy_domains - 1
         ptr_Index_Intrfc(nbvi) = k
         ideb                   = 1
         do i = 1,size_y
           l                    = ideb + (i-1)*size_x
           Index_Intrfc(k)      = l
           k                    = k + 1
           iwork2(l)            = me - xy_domains - 1
         enddo
       endif  
*
       if ((down.eq.1).and.(right.eq.1).and.(bottom.eq.1)) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me - xy_domains - x_domains + 1
         ptr_Index_Intrfc(nbvi) = k
         l                      = size_x
         Index_Intrfc(k)        = l
         k                      = k + 1
         iwork2(l)              = me - xy_domains - x_domains + 1
       endif  
*
       if ((down.eq.1).and.(bottom.eq.1)) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me - xy_domains - x_domains
         ptr_Index_Intrfc(nbvi) = k
         ideb                   = 1
         do i = 1,size_x
           l                    = ideb + (i-1)
           Index_Intrfc(k)      = l
           k                    = k + 1
           iwork2(l)            = me - xy_domains - x_domains
         enddo
       endif   
*
       if ((down.eq.1).and.(left.eq.1).and.(bottom.eq.1)) then
         nbvi                   = nbvi + 1
         indexVi(nbvi)          = me - xy_domains - x_domains - 1
         ptr_Index_Intrfc(nbvi) = k
         l                      = 1
         Index_Intrfc(k)        = l
         k                      = k + 1
         iwork2(l)              = me - xy_domains - x_domains - 1
       endif           
*
*******
       ptr_Index_Intrfc(nbvi+1) = k
*
* Initiate the list of the indices of the interface points
*

c$$$       sizeIntrf = 0
c$$$       do i=1,size_x*size_y*size_z
c$$$         if (iwork(i).ne.0) then
c$$$           sizeIntrf              = sizeIntrf + 1
c$$$           list_Intrfc(sizeIntrf) = i
c$$$           perm(sizeIntrf)        = i
c$$$           weight(sizeIntrf)      = 1.0d0/dble(2**iwork(i))
c$$$         endif
c$$$       enddo
c$$$*
c$$$       j = sizeIntrf+1
c$$$       do i = 1,size_x*size_y*size_z  
c$$$         if (iwork(i).eq.0) then
c$$$           perm(j) = i
c$$$           j = j+1
c$$$         endif
c$$$       enddo


       sizeIntrf = 0
       do i=1,size_x*size_y*size_z
          if (iwork(i).ne.0) then
             sizeIntrf              = sizeIntrf + 1
          endif
       enddo
       
*     
       j = 1
       do i = 1,size_x*size_y*size_z  
          if (iwork(i).eq.0) then
             perm(j) = i
             j = j+1
          endif
       enddo
       
       k=1
       do i=1,size_x*size_y*size_z
          if (iwork(i).ne.0) then
             perm(j)        = i
             j=j+1

             list_Intrfc(k) = i
             weight(k)      = 1.0d0/dble(2**iwork(i))
             k=k+1
          endif
       enddo
       

      Write(me+20,*) "perm_before(:)=",Perm(1:size_x*size_y*size_z)

*
*yoh> 
*     At this point, the schur will be on the "first" rows/columns of the matrix.
*     New maphys version request that the Schur is on the "last" rows/columns of the matrix.
c$$$       Do i=1,sizeIntrf
c$$$          tmp = perm(size_x*size_y*size_z - sizeIntrf +i) 
c$$$          perm(size_x*size_y*size_z - sizeIntrf +i) = perm(i)
c$$$          perm(i) = tmp
c$$$       End Dox
*yoh<


       do i = 1,size_x*size_y*size_z
         invPerm(perm(i)) = i
       enddo
*
* Compute the list of interface points logically owned by the processor
*
       SizeInterf_lo = 0
       do i = 1,size_x*size_y*size_z
         if ((iwork2(i).ne.-1).and.(iwork2(i).gt.me)) then
            SizeInterf_lo = SizeInterf_lo + 1
            list_Intrfc_lo(SizeInterf_lo) = i
         endif
       enddo 
*
       DEALLOCATE(iwork2)
**********************************************************************
* write the local iwork of all entries on a file
**********************************************************************
*      if (me.eq.0) then
*            open(unit=30,file="IWORK",status="REPLACE",iostat=ios)
*            write(unit=30,'(A)'),"%%MatrixMarket matrix coordinate real
*     &      general" 
*            write(unit=30,'(I6,I6,I6)'),size_x*size_y*size_z,1,
*     &      size_x*size_y*size_z
*            write(unit=30,'(I10,I10,I10)'),
*     &      (i,1,iwork(i),i=1,size_x*size_y*size_z)
*            close(unit=30)
*      endif 
***********************************************************************
       return
       end
*
