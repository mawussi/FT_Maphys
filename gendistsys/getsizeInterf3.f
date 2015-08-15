C
C> Get the size of the interface of a domain in a cube [0,1]x[0,1]x[0,1]
C! 
C! @param [in] me          MPI rank of the domain
C! @param [in] x_domains   number of domains in the x direction
C! @param [in] y_domains   number of domains in the y direction
C! @param [in] z_domains   number of domains in the z direction
C! @param [in] size_x      in the domain, number of elements in the x direction
C! @param [in] size_y      in the domain, number of elements in the y direction
C! @param [in] size_z      in the domain, number of elements in the z direction
C! @param [out] iwork      integer workspace of size size_x*size_y*size_z.
C!
C
       integer function getsizeInterf(me,x_domains,y_domains,z_domains,
     &          size_x, size_y,size_z,iwork)
*
* input parameters
       integer me,x_domains,y_domains,z_domains,size_x,size_y,size_z
       integer iwork(*)
* local variables
       integer idom, idom1, sizeIntrf, size_xy, xy_domains
       integer i, j ,k, ideb, l, right, left, top, bottom, up, down, icr
*
*
* setup the default value of weight (ie one everywhere)
*
       idom = me+1
*
       do i=1,size_x*size_y*size_z
          iwork(i) = 0
       enddo
* define the neighboring domains, and sizes
*
       size_xy    = size_x*size_y
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
           right             = 1
           ideb              = size_x
           do j=1,size_z
             do i=1,size_y
               l             = ideb + (i-1)*size_x + (j-1)*size_xy
               iwork(l)      = iwork(l) + 1
               k = k+1
             enddo
           enddo
         endif
         if (mod(idom,x_domains).ne.1) then
* check if the domain has a left neighbor
           left              = 1
           ideb              = 1
           do j=1,size_z
             do i=1,size_y
               l             = ideb + (i-1)*size_x + (j-1)*size_xy
               iwork(l)      = iwork(l) + 1
               k             = k+1
             enddo
           enddo
         endif
       endif
*
       if (y_domains.gt.1) then
* check if the domain has a top neighbor
         idom1 = mod(idom-1,xy_domains)
         if ((idom1/x_domains).ne.(y_domains-1)) then
           top               = 1
           ideb              = (size_y-1)*size_x + 1
           do j=1,size_z
             do i=1,size_x
               l             = ideb + (j-1)*size_xy + (i-1)
               iwork(l)      = iwork(l) + 1
               k             = k+1
             enddo
           enddo
         endif
* check if the domain has a bottom neighbor
         if ((idom1/x_domains).ne.0) then
           bottom            = 1
           ideb              = 1
           do j=1,size_z
             do i=1,size_x
               l             = ideb + (j-1)*size_xy + (i-1)
               iwork(l)      = iwork(l) + 1
               k             = k+1
             enddo
           enddo
         endif
       endif
*
       if (z_domains.gt.1) then
* check if the domain has an up (above) neighbor
         if (((idom-1)/xy_domains).ne.(z_domains-1)) then
           up                = 1
           ideb              = size_xy*(size_z-1)
           do j=1,size_y
             do i=1,size_x
               l             = ideb + (j-1)*size_x + i
               iwork(l)      = iwork(l) + 1
               k             = k+1
             enddo
           enddo
         endif
* check if the domain has a down (below) neighbor
         if (((idom-1)/xy_domains).ne.0) then
           down              = 1
           ideb              = 0
           do j=1,size_y
             do i=1,size_x
               l             = ideb + (j-1)*size_x + i
               iwork(l)      = iwork(l) + 1
               k             = k+1
             enddo
           enddo 
         endif
       endif
*
* Initiate the list of the index of the interface point
*
       sizeIntrf = 0
       do i=1,size_x*size_y*size_z
         if (iwork(i).ne.0) then
           sizeIntrf              = sizeIntrf + 1
         endif
       enddo
*
       getsizeInterf = sizeIntrf
*
       return
       end
