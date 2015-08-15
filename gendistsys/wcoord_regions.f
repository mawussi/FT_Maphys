*     regions problem
*     ac,af,ad are parameters that specify the w_regionsues of respectively
*     c,f,d in figure 4 problem D-F2of 'Local preconditionners for two-level
*     non-overlapping domain decomposition methods' [Carw_regionsho,Giraud,Meurant]
*
      real*8 function w_regions(x,y,z,h1x,h1y,h1z,ac,ad,af)
*     input parameters
      real*4 x,y,z
      real*4 xc,yc,zc
      real*8 ac,ad,af,height,heightOx
      integer h1x,h1y,h1z

      integer comm,me,infompi
      include 'mpif.h'
      comm = MPI_COMM_WORLD
      call MPI_COMM_RANK(comm,me,infompi)

*         if(me.eq.0)write(unit=6,FMT='(F14.3,F14.3,F14.3,F14.3)') 
*     &              x,y,z

*      print*,x,y,z 
*
      height   = 0.2d0
      heightOx = height/2.0d0
      xc = x-1.0
      yc = y-1.0
      zc = z-1.0

*      if(me.eq.0)write(unit=6,
*     &        FMT='(F14.3,F14.3,F14.3,F14.3,F14.3,F14.3,F14.3,F14.3)') 
*     &              x,y,z,xc,yc,zc
*      if(me.eq.0) print*,xc,yc,zc

      if ( (xc.gt.(dble(h1x)*0.4d0)).and.(xc.lt.(dble(h1x)*0.6d0))
     &                .and.(yc.gt.(dble(h1y)*(1.0d0 - heightOx))) ) then
!           oxyde omega1
            w_regions = ac
      elseif ( (xc.lt.(dble(h1x)*0.15d0))
     &        .and.(yc.gt.(dble(h1y)*1.0d0 - height)) ) then
!           top left dopped region omega3
            w_regions = af 
      elseif((xc.gt.(dble(h1x)*0.25d0)).and.(xc.lt.(dble(h1x)*0.45d0))
     &        .and.(yc .gt.(dble(h1y)*(1.0d0-height))) ) then
!           middle left dopped region omega4
            w_regions = af 
      elseif((xc.gt.(dble(h1x)*0.55d0)).and.(xc.lt.(dble(h1x)*0.75d0))
     &        .and.(yc .gt.(dble(h1y)*(1.0d0-height))) ) then
!           middle rigth dopped region omega5
            w_regions = af 
      elseif((xc.gt.(dble(h1x)*0.85d0))
     &        .and.(yc.gt.(dble(h1y)*(1.0d0 - height))) ) then
!           top rigth dopped region omega6
            w_regions = af 
      else !omega 2
            w_regions = ad 
      endif
          end


