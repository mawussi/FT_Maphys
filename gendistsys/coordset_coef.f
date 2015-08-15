C> Generate the coefficients of the matrix 
C!
C! Generate the coefficient of the matrix according to the "problem",
C! and its position
C!
C! @param [in] problem
C!   Indicates the type of system. valid values are :
C!     0 : constant problem
C!     1 : abrupt problem (TODO)
C!     2 : difficult problem (TODO)
C!     3 : sines problem (TODO)
C!     4 : regions problem (TODO)
C!     5 : hd problem with coord hd ?
C!     6-8 : (TODO)
C! @param [in] p1, p2 
C!   Parameters for the convective part.
C! @param [in] an     
C!   Indicates if the system is anisotropic or not.
C!     0 : system is isotropic (ie not anysotropic)
C!     Other : system is anysotorpic
C! @param [in] ac, ad, af 
C!   Parameters for the anisotropic cases.
C!
       subroutine coordset_coef(x,y,z,h1x,h1y,h1z,hx2,hy2,hz2,
     &          s,xm,xp,ym,yp,zm,zp,p1,p2,problem,an,ac,ad,af,
     &                                               hx,hy,hz)


*
*  input parameters
      implicit none
      real*8  hx,hy,hz
      real*4  x,y,z,hx2,hy2,hz2
      real*8  ac,ad,af   ! anisotropic params an==0 <=> no anisotropy
      integer h1x,h1y,h1z
      integer an
      integer problem
*  external functions
      external  w_hd,w_constant,w_regions
      real*8    w_hd,w_constant,w_regions
*  ouput paramaters
      real*8  s, p1,p2
      real*8  xm,xp,ym,yp,zm,zp
*       zp yp
*        |/ 
*    xm--s--xp      zp= z+ coeff, zm = z- coeff
*       /|
*     ym zm

*      ###############################################
*             discretize:  -delta(u) = f
*      ###############################################

*      integer comm,me,infompi
*      include 'mpif.h'
*      comm = MPI_COMM_WORLD
*      call MPI_COMM_RANK(comm,me,infompi)

*         if(me.eq.0)print*,'###########################################'
*         if(me.eq.0)write(unit=6,
*     &        FMT='(F14.3,F14.3,F14.3,F14.3,F14.3)') 
*     &              x,y,z
*         if(me.eq.0)print*,'###########################################'



      p1  = 0
      p2  = 0
*  force -1 on xm and xp for x oriented anisotropy
      if (an.ne.0) then                  ! isotropic
         xm   = -1.0d0
         xp   = -1.0d0
         if (problem.eq.0) then
*           constant problem
            ym   = -w_constant(x   ,y-hy2,z)
            yp   = -w_constant(x   ,y+hy2,z)
            zm   = -w_constant(x   ,y   ,z-hz2)
            zp   = -w_constant(x   ,y   ,z+hz2)
         elseif (problem.eq.1) then
*	    abrupt problem
         elseif (problem.eq.2) then
*	    difficult problem
         elseif (problem.eq.3) then
*	    sines problem
         elseif (problem.eq.4) then
*	    regions problem
            ym   = -w_regions(x   ,y-hy2,z   ,h1x,h1y,h1z,ac,ad,af)
            yp   = -w_regions(x   ,y+hy2,z   ,h1x,h1y,h1z,ac,ad,af)
            zm   = -w_regions(x   ,y   ,z-hz2,h1x,h1y,h1z,ac,ad,af)
            zp   = -w_regions(x   ,y   ,z+hz2,h1x,h1y,h1z,ac,ad,af)

         elseif (problem.eq.5) then
*	    hd problemwcoord_hd(x,y,z,h1x,h1y,h1z)
            ym   = -w_hd(x   ,y-hy2,z   ,h1x,h1y,h1z)
            yp   = -w_hd(x   ,y+hy2,z   ,h1x,h1y,h1z)
            zm   = -w_hd(x   ,y   ,z-hz2,h1x,h1y,h1z)
            zp   = -w_hd(x   ,y   ,z+hz2,h1x,h1y,h1z)
         elseif (problem.eq.6) then
*           porous media
         elseif (problem.eq.7) then
*           pyramid media
         elseif (problem.eq.8) then
*           3D regions semi conductor
         endif
      else                        ! anisotropic ?
         if (problem.eq.0) then
*           constant problem
            xm   = -w_constant(x-hx2,y   ,z)
            xp   = -w_constant(x+hx2,y   ,z)
            ym   = -w_constant(x   ,y-hy2,z)
            yp   = -w_constant(x   ,y+hy2,z)
            zm   = -w_constant(x   ,y   ,z-hz2)
            zp   = -w_constant(x   ,y   ,z+hz2)
         elseif (problem.eq.1) then
*	    abrupt problem
         elseif (problem.eq.2) then
*	    difficult problem
         elseif (problem.eq.3) then
*	    sines problem
         elseif (problem.eq.4) then
*           regions problem

*         if(me.eq.0)write(unit=6,
*     &        FMT='(F14.3,F14.3,F14.3,F14.3,F14.3,F14.3,F14.3,F14.3)') 
*     &              x,y,z,x-hx2,y,z
            xm   = -w_regions(x-hx2,y   ,z,   h1x,h1y,h1z,ac,ad,af)
*         if(me.eq.0)write(unit=6,
*     &        FMT='(F14.3,F14.3,F14.3,F14.3,F14.3,F14.3,F14.3,F14.3)') 
*     &              x,y,z,x+hx2,y,z
            xp   = -w_regions(x+hx2,y   ,z,   h1x,h1y,h1z,ac,ad,af)
*         if(me.eq.0)write(unit=6,
*     &        FMT='(F14.3,F14.3,F14.3,F14.3,F14.3,F14.3,F14.3,F14.3)') 
*     &              x,y,z,x,y-hy2,z
            ym   = -w_regions(x   ,y-hy2,z,   h1x,h1y,h1z,ac,ad,af)
*         if(me.eq.0)write(unit=6,
*     &        FMT='(F14.3,F14.3,F14.3,F14.3,F14.3,F14.3,F14.3,F14.3)') 
*     &              x,y,z,x,y+hy2,z
            yp   = -w_regions(x   ,y+hy2,z,   h1x,h1y,h1z,ac,ad,af)
*         if(me.eq.0)write(unit=6,
*     &        FMT='(F14.3,F14.3,F14.3,F14.3,F14.3,F14.3,F14.3,F14.3)') 
*     &              x,y,z,x,y,z-hz2
            zm   = -w_regions(x   ,y   ,z-hz2,h1x,h1y,h1z,ac,ad,af)
*         if(me.eq.0)write(unit=6,
*     &        FMT='(F14.3,F14.3,F14.3,F14.3,F14.3,F14.3,F14.3,F14.3)') 
*     &              x,y,z,x,y,z+hz2
            zp   = -w_regions(x   ,y   ,z+hz2,h1x,h1y,h1z,ac,ad,af)
*         call sleep(1)
         elseif (problem.eq.5) then
*	    hd problem
            xm   = -w_hd(x-hx2,y    ,z    ,h1x,h1y,h1z)
            xp   = -w_hd(x+hx2,y    ,z    ,h1x,h1y,h1z)
            ym   = -w_hd(x    ,y-hy2,z    ,h1x,h1y,h1z)
            yp   = -w_hd(x    ,y+hy2,z    ,h1x,h1y,h1z)
            zm   = -w_hd(x    ,y    ,z-hz2,h1x,h1y,h1z)
            zp   = -w_hd(x    ,y    ,z+hz2,h1x,h1y,h1z)
         elseif (problem.eq.6) then
*           porous media
         elseif (problem.eq.7) then
*           pyramid media
         elseif (problem.eq.8) then
*           3D regions semi conductor 
         endif
      endif
      s   = -xm -ym -xp -yp - zm - zp
      return
      end
