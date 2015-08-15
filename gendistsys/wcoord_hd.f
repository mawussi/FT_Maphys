      real*8 function w_hd(x,y,z,h1x,h1y,h1z)
*     input parameters
      integer h1x,h1y,h1z
      real*4  x,y,z
      real*8  xc1,xc2,yc1

      integer comm,me,nproc,infompi
      include 'mpif.h'
      comm = MPI_COMM_WORLD
      call MPI_COMM_SIZE(comm,nproc,infompi)
      call MPI_COMM_RANK(comm,me,infompi)



      xc1= dble(h1x)/3.0d0
      xc2= dble((h1x)*2.0d0)/3.0d0
      yc1= dble(h1y)/2.0d0

*	return value 
        if ((y-1.0) .lt. yc1) then
          if ( ((x-1.0) .lt. xc1) .or. ((x-1.0) .gt. xc2) ) then
             w_hd = 1.0d0
          else
             w_hd = 1.0d3
          endif
        else
          if ( ((x-1.0) .lt. xc1) .or. ((x-1.0) .gt. xc2) ) then
             w_hd = 1.0d3
          else
             w_hd = 1.0d0
          endif
        endif




*         if(me.eq.0)write(unit=6,
*     &      FMT='(F14.3,F14.3,F14.3,F14.3,F14.3,F14.3,F14.3)') 
*     &             xc1,xc2,yc1, x,y,z,wcoord_hd

        return
        end

