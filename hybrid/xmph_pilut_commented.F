c Warning: XMPH_GENFILE_COMMENT
c----------------------------------------------------------------------
c
c Partial ILUT - with Azzam modifications 
c
c----------------------------------------------------------------------
c Changelog :
c  2011-10-06 : Yohan Lee-tin-yien 
c               Remove all useless arguments
c               Add the statistics to be printed if STAT_LVL is defined
c  2011-09-15 : Yohan Lee-tin-yien plug into new MaPHyS
c  before 2011-09-15 : Azzam Modifications to MaPHyS prototype
c  before Azzam modfications : Initial implementation 
c
c----------------------------------------------------------------------
c cori is the code in the original version, which we commented.
c Clg  are the comments added by Luc Giraud
c Cyl  are the comments added by Yohan Lee-tin-yien
c Cah  are the comments added by Azzam Haidar
c----------------------------------------------------------------------
c Preprocessing
c----------------------------------------------------------------------

Cyl define MPH_INT
#include "mph_defs_f.h"
#include "mph_macros_f.h"

Cyl CONVERT2FLOAT
#if XMPH_HAVE_ARITH_C
#define CONVERT2FLOAT(x) CMPLX(x,KIND=4)       
#elif XMPH_HAVE_ARITH_D
#define CONVERT2FLOAT(x) REAL(x,KIND=8)       
#elif XMPH_HAVE_ARITH_S
#define CONVERT2FLOAT(x) REAL(x,KIND=4)       
#elif XMPH_HAVE_ARITH_Z
#define CONVERT2FLOAT(x) CMPLX(x,KIND=8)       
#endif

Cyl controls the level of statistics of the routine
#define STAT_PRINT 1
#define STAT_EXPORT 2

#ifndef STAT_LVL
#define STAT_LVL 0
#endif

#if STAT_LVL >= STAT_EXPORT
#define IF_STAT_EXPORT(x) x
#else
#define IF_STAT_EXPORT(x) !x
#endif

c----------------------------------------------------------------------- 
c PARTIAL ILUT - adapted from ILUT. Mon Aug 24 11:48:10 CDT 1998 YS. 
c-----------------------------------------------------------------------
cori 
cori      SUBROUTINE XMPH_pilut
cori     *    (n,a,ja,ia,lfil,droptol,last,iptr,lalu,alu,jlu,ju,
cori     *     lxm,xm,jxm,ixm,w,jw,ierr)

cyl      SUBROUTINE XMPH_pilut
cyl     *    (n,a,ja,ia,lfil,droptol,last,iptr,lalu,alu,jlu,ju,
cyl     *     w,jw,ierr,lfilschur,dropschur,
cyl     &     fullU_entries,MYRINFO)

      SUBROUTINE XMPH_pilut
     *     (n,a,ja,ia,lfil,droptol,last,lalu,alu,jlu,ju,
     *     w,jw,ierr,lfilschur,dropschur,tpilu)

#if STAT_LVL >= STAT_PRINT
      Use mph_dbg_mod
#endif
      IMPLICIT NONE 

*     input arguments
      MPH_INT n,lfil,last,lalu
      MPH_INT lfilschur
      REAL*8 droptol,dropschur
      MPH_INT ja(*),ia(n+1)
      XMPH_FLOAT a(*)


*     output arguments
      INTEGER ierr
      MPH_INT ju(2*n+1)
      MPH_INT jlu(*) 
      XMPH_FLOAT alu(*)
      Real*8 tpilu

*     in/out arguments
      MPH_INT jw(2*n)
      XMPH_FLOAT w(n)

*     obsolete & to suppress arguments
cori  MPH_INT iptr
cori  MPH_INT lxm
cori  MPH_INT jxm(*),ixm(*) 
cori  XMPH_FLOAT xm(*)
cyl   REAL*8 MYRINFO(*)
cyl   MPH_INT fullU_entries

c----------------------------------------------------------------------
c --> this version does not drop outside the schur complement.
c----------------------------------------------------------------------- 
c Partial block ILU factorization with dual truncation. 
c Consider a matrix partitioned and factored as follows 
c                                                                      
c |  B  F |        | L    0 |   | U   U1 |
c |  E  C |   =    | L1   I | * | 0   S  |                   
c                                                                      
c where B is a sub-matrix of dimension=last
c 
c piliut computes the  above *partial* imcomplete factorization.
c LU ~= B, L1 ~= E U^{-1}, etc..  
c From 1 to last: it  will proceed as in ILUT
c From rows  last to n :  it will  process  the rows as in  ILUT --
c except  that the columns  elements elminated are  only those from 1 to
c last. The result is   that  the (2,2)  matrix   S obtained from   this
c factorization  is  exactly the Schur  complement.    This code uses an
c adapted version of  ILUT (which is based on  CSR format instead of the
c MSR format.) This is coded to take advantage of the existing structure
c of ILUT.
C
Clg The initial routine is by Y. Saad for Arms. We change it to compute
Clg an sparse approximation of the Schur complement based on an ILU 
Clg factorization of the (1,1) block B.
Clg The ILU factorization is based on a right-looking approach. However,
Clg because the matrix is given in a CSR format, it should be called
Clg a "up-looking" variant.
Clg As a consequence, because we are only interested in the approximated
Clg Schur complement, the entries computed in the L part are not kept.
Clg In addition, we used two different thresholds, the first "droptol" is 
Clg the regular ILU threshold in the U and U1 matrices. The second one is 
Clg "dropschur", that is used for the droppting in the approximated Schur.
Clg Similarly and extra fil parameter has been added to enable more fill
Clg in the Schur complement than in the U and U1 parts.
c
c----------------------------------------------------------------------
c
c on entry:
c========== 
c n       = The row dimension of the matrix A. 
c
c a,ja,ia = input matrix stored in Compressed Sparse Row format.              
c
c lfil    = The fill-in parameter. Same meaning as in ILUT. On return,
c           each row of the L-part and each row of the U-part will have 
c           lfil nonzero entries at most (excluding the diagonal 
c           elements). Note that the L-part and U-part are meant in the
c           above block sense. [the U part consists of the rows of the
c           matrices U, U1  and S, The L-part consists of the rows of
c           L and L1 ].    lfil must be .ge. 0.
c
c droptol =  Sets the threshold for dropping small terms in the
c           factorization. See ILUT for details. 
c
c last    = Number of the last row to be fully processed. Equal to
c           the dimension of the  block B in above illustration. 
c
c iptr    = integer. Starting point of the current level in (xm,jxm,ixm) 
c  
c lalu    = integer. The lengths of arrays alu and jlu -- see below.
c           lalu should be at least: (2*lfil + 1)*n
c           Note that ju should also be of length 2*n+1. 
c
cori lxm     = length of arrays xm and jxm that is left. Remember 
cori           that the useful part of xm, jxm starts in location iptr.
cori           lxm refers to the total length of xm, jxm. The space that is
cori           consumed in pilut is :  (n + last)*lfil + last.
cori           (so lxm should be .ge. (iptr-1) + (n + last)*lfil + last 
c
c
c lfilschur 
c         = The fill-in parameter of the schur complement. On return,
c           each row of the schur complement will have 
c           lfil nonzero entries at most (excluding the diagonal 
c           elements).
c           (variables added by Azzam, comments added by yohan)
c
c schurtol=  Sets the threshold for dropping small terms in the
c           schur complement.
c           (added by Azzam, commented by yohan) 
c 
c fullU_entries 
c         =  Say how many entries there are in a full LU factorisation.
c            This is used by the routine to print several 
c            statistics.
c           (variables added by Azzam, comments added by yohan)
c             
c On return:
c===========
c
c alu,jlu,
c ju      = contains the Schur complement matrix (S in above diagram)
c           stored in CSR format. 
c          NOTE: ju is initially used as a work array and its length
c          is assumed to be 2*n+1 on entry.
c
cori xm, jxm, ixm, 
cori         = contains (in sequence) two matrices both stored in CSR
cori           format -- one after the other: 1) the matrix [U, U1] 
cori           of size last above diagram. 2) The matrix of size n 
cori           consisting of  | L  |  of the above diagram (ones not 
cori                          | L1 |  represented.) 
c
c ierr    = integer. Error message with the following meaning.
c           ierr  =  0   --> successful return.
c           ierr  =  1   --> Error. input matrix may be wrong.
c                            (The elimination process has generated a
c                            row in L or U whose length is .gt.  n.)
c           ierr  =  2   --> The matrix L overflows the arrays alu,jlu 
c           ierr  =  3   --> The matrix U overflows the arrays alu,jlu 
c           ierr  =  4   --> The array xm overflows 
c           ierr  =  5   --> Illegal value for lfil or last 
c           ierr  =  6   --> zero row encountered.
c           ierr  =  7   --> Tnorm divided by zero
c                            (Azzam Haidar, Yohan Lee-tin-yien)
c
c----------------------------------------------------------------------- 
c work arrays:
c=============
c jw        = integer work array of length 2*n
c w         = real work array of length n. 
c ju        = is also used as work array. Length should be 2*n+1.  
c  
c----------------------------------------------------------------------
c w, ju(1:n) store the working array [1:ii-1 = L-part, ii:n = u] 
c     locals
      MPH_INT ju0,k,j1,j2,j,ii,i,lenl,lenu,jj,jrow,jpos,len, jshift 
      MPH_INT next, ibound, jcol, nsch, ier
cori  MPH_INT idxm
      REAL*8  tnorm
cyl   REAL*8  t ! previously commented in lines 574,561. Seems useless
      XMPH_FLOAT tempp,s,fact

c----------------------------------------------------------------------- 
c     variables added from original version
c----------------------------------------------------------------------- 

#if STAT_LVL >= STAT_PRINT
      MPH_INT nnzU,nnzS,nnzUS,lenuinit
      MPH_INT memoU,memoUinit

      REAL*8 keepU,keeplfilU,keepS,keeplfilS
      REAL*8 minpercentU,minpercentlfilU 
      MPH_INT ligneU,lignelfilU
      MPH_INT size_S
      MPH_INT keptdropvalU,keptlfilvalU,fullvalU
      MPH_INT keptdropvalS,keptlfilvalS,fullvalS
#endif
#if STAT_LVL >= STAT_EXPORT
      REAL*8, dimension(:,:),allocatable :: percentage
#endif 
*     constants 
      Real(Kind=8), Parameter :: DZERO = 0._8
*     timing
      Real(Kind=8), External :: MPI_Wtime


c-End of header---------------------------------------------------------


      if (lfil .lt. 0 .or. last .le. 0) goto 998
c-----------------------------------------------------------------------
c     initialize ju0 (points to next element to be added to alu,jlu)
c     and pointer array.
c-----------------------------------------------------------------------
      ju0 = 1
      ju(n+1) = ju0 
c
c     initialize nonzero indicator array. 
c
      do 1 j=1,n
         jw(n+j)  = 0
 1    continue
c
c     beginning of main loop.
c
#if STAT_LVL >= STAT_PRINT
      Call MPH_dbg_init()
      size_S=0
      fullvalU     = 0
      keptdropvalU = 0
      keptlfilvalU = 0
      fullvalS     = 0
      keptdropvalS = 0
      keptlfilvalS = 0

      keepU=0.0d0
      keepS=0.0d0
      keeplfilU=0.0d0
      keeplfilS=0.0d0

      ligneU=0
      lignelfilU=0
      minpercentU=1.0d3
      minpercentlfilU=1.0d3

      memoU     = lfil*last 
      memoUinit = memoU 
#endif
#if STAT_LVL >= STAT_EXPORT      
      ALLOCATE(percentage(2,n))
      percentage(:,:)=0.0d0
#endif
      tpilu = MPI_Wtime()
      next = last+1
      jshift = 0 
      do 500 ii = 1, n
         j1 = ia(ii)
         j2 = ia(ii+1) - 1
         tnorm = DZERO
         do 501 k=j1,j2
            tnorm = tnorm+ABS(a(k))
 501     continue
         if (tnorm .eq. DZERO)then
            if(ii.le.last) goto 999
            ! ii>last this is a Schur line 
            ibound = ii
            len=0
            goto 331  
         endif

cyoh>     check added by Azzam Haidar
         if((j2-j1+1).eq.0) Goto 993
cyoh<
         tnorm = tnorm/dble(j2-j1+1)
         ibound = min(ii,next)
         if (ii .gt. last) jshift = last
c     
c     unpack L-part and U-part of row of A in arrays w 
c     
         lenu = 1
         lenl = 0
c
         jw(ibound) = ii
         w(ibound) = XMPH_FLOATZERO
         jw(n+ii) = ibound 
c
         do 170  j = j1, j2
            jcol = ja(j)
            tempp = a(j)
            if (jcol .lt. ibound) then
               lenl = lenl+1
               jw(lenl) = jcol
               w(lenl) = tempp
               jw(n+jcol) = lenl
            else if (jcol .eq. ii) then
               w(ibound) = tempp
            else
               lenu = lenu+1
               jpos = ibound+lenu-1 
               jw(jpos) = jcol
               w(jpos) = tempp
               jw(n+jcol) = jpos
            endif
 170     continue
         jj = 0
         len = 0 
c     
c     eliminate previous rows
c     
 150     jj = jj+1
         if (jj .gt. lenl) goto 160
c-----------------------------------------------------------------------
c     in order to do the elimination in the correct order we must select
c     the smallest column index among jw(k), k=jj+1, ..., lenl.
c-----------------------------------------------------------------------
         jrow = jw(jj)
         k = jj
c     
c     determine smallest column index
c     
         do 151 j=jj+1,lenl
            if (jw(j) .lt. jrow) then
               jrow = jw(j)
               k = j
            endif
 151     continue
c
         if (k .ne. jj) then
c     exchange in jw
            j = jw(jj)
            jw(jj) = jw(k)
            jw(k) = j
c     exchange in jr
            jw(n+jrow) = jj
            jw(n+j) = k
c     exchange in w
            s = w(jj)
            w(jj) = w(k)
            w(k) = s
         endif
c
c     zero out element in row by setting jw(n+jrow) to zero.
c     
         jw(n+jrow) = 0
c
c     get the multiplier for row to be eliminated (jrow).
c     
         fact = w(jj)*alu(ju(jrow))
cyoh  In original version, small pivots were avoided
cori         if ((abs(fact) .le. droptol).and.(ii.le.last)) goto 150

c     
c     combine current row and row jrow
c
         do 203 k = ju(jrow)+1, ju(n+jrow+1)-1
            s = fact*alu(k)
            j = jlu(k)
            jpos = jw(n+j)
            if (j .ge. ibound) then
c     
c     dealing with upper part.
c     
               if (jpos .eq. 0) then
c     
c     this is a fill-in element
c     
                  lenu = lenu+1
                  if (lenu .gt. n) goto 994 ! instaed of "n" it should be "(n-ii+1)" for ii.le.last 
                                            ! and "(n-last+1)" for ii> last
                  i = ibound +lenu-1
                  jw(i) = j
                  jw(n+j) = i
                  w(i) = - s
               else
c
c     this is not a fill-in element 
c
                  w(jpos) = w(jpos) - s
               endif
            else
c     
c     dealing  with lower part.
c     
               if (jpos .eq. 0) then
c     
c     this is a fill-in element
c     
                  lenl = lenl+1
                  if (lenl .gt. (ii-1)) goto 994
                  jw(lenl) = j
                  jw(n+j) = lenl
                  w(lenl) = - s
               else
c     
c     this is not a fill-in element 
c     
                  w(jpos) = w(jpos) - s
               endif
            endif
 203     continue
c     
c     store this pivot element -- (from left to right -- no danger of
c     overlap with the working elements in L (pivots). 
c     
         len = len+1 
         w(len) = fact
         jw(len)  = jrow
         goto 150
 160     continue
c     
c     reset double-pointer to zero (U-part)
c     
         do 308 k=1, lenu
            jw(n+jw(ibound+k-1)) = 0
 308     continue
c     
c     update L-matrix
c     
cori         lenl = len 
cori         len = min0(lenl,lfil)
c     
c     sort by quick-split
c     
cori         call qsplit (w,jw,lenl,len)
c
c     store L-part
c 
c ?? 
cori         do 204 k=1, len 
cori            if (ju0 .gt. lalu) goto 995
cori            alu(ju0) =  w(k)
cori            jlu(ju0) =  jw(k)
cori            ju0 = ju0+1
cori 204     continue
c
c     update U-matrix -- first apply dropping strategy 
c
         if(ii.le.last) then  ! U part 
            len = 0
            do k=1, lenu-1
               if (ABS(w(ibound+k)) .gt. droptol*tnorm) then 
                  len = len+1
                  w(ibound+len) = w(ibound+k) 
                  jw(ibound+len) = jw(ibound+k) 
               endif
            enddo

#if STAT_LVL >= STAT_PRINT
            fullvalU=fullvalU+lenu 
            if(len.gt.0) then 
               keepU=((dble(len+1))/dble(lenu))*100.0d0
               IF_STAT_EXPORT( percentage(1,ii)=keepU ) 
               minpercentU = min(minpercentU,keepU) 
               if(minpercentU.eq.keepU) ligneU=ii 
               keptdropvalU=keptdropvalU+len+1
            endif
            lenuinit=lenu
#endif
            lenu = len+1
            len = min0(lenu,lfil)

#if STAT_LVL >= STAT_PRINT 
            keeplfilU=(dble(len)/dble(lenuinit))*100.0d0
            IF_STAT_EXPORT(percentage(2,ii)= keeplfilU)
            minpercentlfilU=min(minpercentlfilU,keeplfilU)
            if(minpercentlfilU.eq.keeplfilU) lignelfilU=ii
            keptlfilvalU=keptlfilvalU+len
            memoU  = memoU-len
#endif

         else ! Schur part U
            len = 0
            do k=1, lenu-1
               if (ABS(w(ibound+k)) .gt. dropschur*tnorm) then 
                  len = len+1
                  w(ibound+len) = w(ibound+k) 
                  jw(ibound+len) = jw(ibound+k) 
               endif
            enddo

#if STAT_LVL >= STAT_PRINT 
            fullvalS=fullvalS+lenu
            if(len.gt.0) then
               keepS=((dble(len+1))/dble(lenu))*100.0d0
               IF_STAT_EXPORT( percentage(1,ii)=keepS )
               keptdropvalS=keptdropvalS+len+1
            endif
            lenuinit=lenu
#endif

            lenu = len+1
            len = min0(lenu,lfilschur)

#if STAT_LVL >= STAT_PRINT 
            keeplfilS=(dble(len)/dble(lenuinit))*100.0d0
            IF_STAT_EXPORT( percentage(2,ii)= keeplfilS )
            keptlfilvalS=keptlfilvalS+len
#endif

         endif
c
         call qsplit (w(ibound+1), jw(ibound+1), lenu-1,len)
c
c     save pointer to diagonal element (=beginning of row ii of U)
c     
 331     ju(ii) = ju0
c     
c     store inverse of diagonal element of u
c     
*         t = ABS(w(ibound))
         if (jshift .le. 0) then 
            if (w(ibound) .eq. XMPH_FLOATZERO) 
     &      w(ibound) = CONVERT2FLOAT((0.0001+droptol)*tnorm)
            w(ibound) = XMPH_FLOATONE/w(ibound) 
         endif
c     
c     do the copying + calculate norm
c     
         if (len + ju0 .gt. lalu) goto 996
         do 302 k=ibound,ibound+len-1 
            jlu(ju0) = jw(k) - jshift 
            alu(ju0) = w(k)
*            t = t + ABS(w(k) )
            ju0 = ju0+1
 302     continue
c     
c     update pointer to beginning of next row of U.
c     
         ju(n+ii+1) = ju0
c-----------------------------------------------------------------------
c     end main loop
c-----------------------------------------------------------------------
 500  continue

      tpilu = MPI_Wtime()-tpilu

#if STAT_LVL >= STAT_PRINT
      nnzU  = (ju(n+last+1)-1)
      nnzS  = (ju(2*n+1)-1)-(ju(n+last+1)-1)
      nnzUS = (ju(2*n+1)-1)
cyl> MYRINFO deleted
cyl      MYRINFO(10)=(dble(keptdropvalU)*100.0d0)/dble(fullvalU)
cyl      MYRINFO(11)=(dble(keptlfilvalU)*100.0d0)/dble(fullvalU)
cyl      MYRINFO(12)=(dble(keptdropvalS)*100.0d0)/dble(fullvalS)
cyl      MYRINFO(13)=(dble(keptlfilvalS)*100.0d0)/dble(fullvalS)

      Call MPH_dbg_set_file("pilut_stats.txt")

      write(unit=dbg_unit,FMT='(A)')
     &      '                                               INFO'
      write(unit=dbg_unit,FMT='(A)')
     &'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write(unit=dbg_unit,FMT='(A,F10.2,A4)')
     &   ' Allocated real space for U                :',
     &   (dble(memoUinit)*8.0d0)/1.0d6,' MB '

         write(unit=dbg_unit,FMT='(A,F10.2,A4)')
     &   ' Used      real space for U                :',
     &   (dble(memoUinit-memoU)*8.0d0)/1.0d6,' MB '

         write(unit=dbg_unit,FMT='(A,F10.2,A4)')
     &   ' Kept      real space for U                :',
     &   (dble(memoU)*8.0d0)/1.0d6,' MB '

         write(unit=dbg_unit,FMT='(A,I10)')
     &   ' Optimal lfil value could be               :',
     &   ((memoUinit-memoU)/last)+1

      write(unit=dbg_unit,FMT='(A)')
     &'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write(unit=dbg_unit,FMT='(A5,I5,A45,3A15)')
     &   '---->',dbg_rank,' ','U','Schur','U+Schur'
         write(unit=dbg_unit,FMT='(A5,I5,A45,3I15)')
     &   '---->',dbg_rank,'  NNZ OF U, Schur, U+Schur :',
     &   nnzU,nnzS,nnzUS
         write(unit=dbg_unit,FMT='(A55,F12.2,A3,F12.2,A3,F12.2,A3)')
     &   ' Required Memory (real) for U, Schur, U+Schur :',
     &   (dble(nnzU)*8.0d0)/1.0d6,' MB ',
     &   (dble(nnzS)*8.0d0)/1.0d6,' MB ',
     &   (dble(nnzUS)*8.0d0)/1.0d6,' MB '
         write(unit=dbg_unit,FMT='(A55,F12.2,A3,F12.2,A3,F12.2,A3)')
     &   ' Required Memory (integer) for U, Schur, U+Schur :',
     &   (dble(nnzU)*4.0d0)/1.0d6,' MB ',
     &   (dble(nnzS)*4.0d0)/1.0d6,' MB ',
     &   (dble(nnzUS)*4.0d0)/1.0d6,' MB '
         write(unit=dbg_unit,FMT='(A55,F12.2,A3,F12.2,A3,F12.2,A3)')
     &   ' Required Memory (all internal) for U, Schur, U+Schur :',
     &   (dble(nnzU)*12.0d0)/1.0d6,' MB ',
     &   (dble(nnzS)*12.0d0)/1.0d6,' MB ',
     &   (dble(nnzUS)*12.0d0)/1.0d6,' MB '
         write(unit=dbg_unit,FMT='(A)')
     &'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write(unit=dbg_unit,FMT='(A,F7.2,A2)')
     &   ' Kept after dropping in Schur              :', keepS,'%'
         write(unit=dbg_unit,FMT='(A,F7.2,A2)')
     &   ' Kept after lfil  S in Schur               :', keeplfilS,'%'
      write(unit=dbg_unit,FMT='(A)')
     &'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(unit=dbg_unit,FMT='(A)') " "
      write(unit=dbg_unit,FMT='(A)') " "
      write(unit=dbg_unit,FMT='(A)')
     &'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(unit=dbg_unit,FMT='(A)')
     &'         AZZAM MODIFICATION 13/09/2009 :'
      write(unit=dbg_unit,FMT='(A)')
     &'         GOOD VALUES USE FULL LU entries as input'  
      write(unit=dbg_unit,FMT='(A)')
     &'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         write(unit=dbg_unit,FMT='(A,F7.2,A2,I10)')
     &   ' MIN % after dropping in U factor           :',
     &   minpercentU,'%',ligneU
         write(unit=dbg_unit,FMT='(A,F7.2,A2,I10)')
     &   ' MIN % after lfil  in U factor              :', 
     &   minpercentlfilU,'%',lignelfilU

cyl> those statistics depends on fullU_entries which is an argument 
cyl  that I removed (since it depends on the sparse direct solver).
cyl
cyl         write(unit=dbg_unit,FMT='(A,4I16)')
cyl     &   'ENTRIES: full_LU, full_ILU, afterdrop, afterlfil :',
cyl     &   fullU_entries,fullvalU,keptdropvalU,keptlfilvalU
cyl         write(unit=dbg_unit,FMT='(A49,A7,A10)')
cyl     &   ' ','%full_LU','%full_ILU' 
cyl         write(unit=dbg_unit,FMT='(A,F7.2,A2,F10.2,A2)')
cyl     &   ' Good Assumed full-ILU percentage of full-LU    :', 
cyl     &   (dble(fullvalU)/dble(fullU_entries))*100.0d0,'%',
cyl     &   (dble(fullvalU)/dble(fullvalU))*100.0d0,'%'
cyl         write(unit=dbg_unit,FMT='(A,F7.2,A2,F10.2,A2)')
cyl     &   ' Good Kept after dropping in U factor           :', 
cyl     &   (dble(keptdropvalU)/dble(fullU_entries))*100.0d0,'%',
cyl     &   (dble(keptdropvalU)/dble(fullvalU))*100.0d0,'%'
cyl         write(unit=dbg_unit,FMT='(A,F7.2,A2,F10.2,A2)')
cyl     &   ' Good Kept after lfil  in U factor              :', 
cyl     &   (dble(keptlfilvalU)/dble(fullU_entries))*100.0d0,'%',
cyl     &   (dble(keptlfilvalU)/dble(fullvalU))*100.0d0,'%'
cyl<
         size_S=n-last
         write(unit=dbg_unit,FMT='(A,4I16)')
     &   'ENTRIES: full_LUS, full_ILUS, afterdrop, afterlfil :',
     &   size_S**2,fullvalS,keptdropvalS,keptlfilvalS
         write(unit=dbg_unit,FMT='(A,F7.2,A2,F10.2,A2)')
     &   ' Good Assumed ILU-S percentage of  EXACT-S      :', 
     &   (dble(fullvalS)/dble(size_S**2))*100.0d0,'%',
     &   (dble(fullvalS)/dble(fullvalS))*100.0d0,'%'
         write(unit=dbg_unit,FMT='(A,F7.2,A2,F10.2,A2)')
     &   ' Good Kept after dropping in Schur              :', 
     &   (dble(keptdropvalS)/dble(size_S**2))*100.0d0,'%',
     &   (dble(keptdropvalS)/dble(fullvalS))*100.0d0,'%'
         write(unit=dbg_unit,FMT='(A,F7.2,A2,F10.2,A2)')
     &   ' Good Kept after lfil  S in Schur               :', 
     &   (dble(keptlfilvalS)/dble(size_S**2))*100.0d0,'%',
     &   (dble(keptlfilvalS)/dble(fullvalS))*100.0d0,'%'
         write(unit=dbg_unit,FMT='(A)')
     &'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      write(unit=dbg_unit,FMT='(A)') " "
#if STAT_LVL >= STAT_EXPORT
      Call MPH_dbg_set_file("pilut_v0_percentages.txt")
Cah      write(unit=dbg_unit,FMT='(2F15.4)') dble(-10), dble(n)
Cah      write(unit=dbg_unit,FMT='(2F15.4)') dble(-10), dble(last)
Cah      write(unit=dbg_unit,FMT='(2F15.4)') dble(-10), dble(n-last)
      write(unit=dbg_unit,FMT='(2F15.4)') 
     &     (percentage(1,i),percentage(2,i),i=1,n) 
      Close(unit=dbg_unit)
#endif
      Call MPH_dbg_exit
#endif

c
c     post-process to get desired data structure --
c     -   copy U-part to beginning of xm array 
cori   idxm = iptr 
cori   call pcopmat(last,alu,jlu,ju,ju(n+2),xm,jxm,ixm,idxm,1,lxm,ier)
cori   if (ier .ne. 0) goto 997
c
c     - copy L-part after 
cori   nsch = n - last
cori   idxm = ixm(next)
cori   call pcopmat(n,alu,jlu,ju(n+1),ju,xm,jxm,ixm(next),idxm,1,
cori  *    lxm,ier)
cori   if (ier .ne. 0) goto 997
c
c     - copy all schur part in alu-jlu array - 
      nsch = n - last
      call pcopmat(nsch,alu,jlu,ju(next),ju(n+next+1),alu,jlu,jw,1,1,
     *     lalu,ier) 
c
c     copy pointer array from jw for schur-part 
c
      do j=1, nsch+1 
         ju(j) = jw(j) 
      enddo
c
c-----done--------------------------------------------------------------
c     correct return. 
      ierr = 0
      return
c     tnorm divided by zero (added by Azzam Haidar)
 993  ierr = 7
      return
c     incomprehensible error. Matrix must be wrong.
 994  ierr =  1
      return
c     insufficient storage in L.
cori 995  ierr =  2
cori      return
c     insufficient storage in U.
 996  ierr =  3
      return
c     array xm overflows. 
cori 997  ierr =  4
cori      return 
c     illegal value for lfil or last entered.
 998  ierr =  5
      return
c     zero row encountered
 999  ierr =  6 
      return
c----------------end-of-ilutp-------------------------------------------
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
        subroutine qsplit(a,ind,n,ncut)
        XMPH_FLOAT a(n)
        MPH_INT ind(n), n, ncut
c-----------------------------------------------------------------------
c     does a quick-sort split of a real array.
c     on input a(1:n). is a real array
c     on output a(1:n) is permuted such that its elements satisfy:
c
c     abs(a(i)) .ge. abs(a(ncut)) for i .lt. ncut and
c     abs(a(i)) .le. abs(a(ncut)) for i .gt. ncut
c
c     ind(1:n) is an integer array which permuted in the same way as a(*).
c-----------------------------------------------------------------------
        XMPH_FLOAT tmp 
        REAL*8 abskey
        MPH_INT itmp, first, last
c-----
        first = 1
        last = n
        if (ncut .lt. first .or. ncut .gt. last) return
c
c     outer loop -- while mid .ne. ncut do
c
 1      mid = first
        abskey = ABS(a(mid))
        do 2 j=first+1, last
           if (ABS(a(j)) .gt. abskey) then
              mid = mid+1
c     interchange
              tmp = a(mid)
              itmp = ind(mid)
              a(mid) = a(j)
              ind(mid) = ind(j)
              a(j)  = tmp
              ind(j) = itmp
           endif
 2      continue
c
c     interchange
c
        tmp = a(mid)
        a(mid) = a(first)
        a(first)  = tmp
c
        itmp = ind(mid)
        ind(mid) = ind(first)
        ind(first) = itmp
c
c     test for while loop
c
        if (mid .eq. ncut) return
        if (mid .gt. ncut) then
           last = mid-1
        else
           first = mid+1
        endif
        goto 1
c----------------end-of-qsplit------------------------------------------
c-----------------------------------------------------------------------
        end
c-----------------------------------------------------------------------
c----------------------------------------------------------------------- 
      subroutine pcopmat (nrow,a,ja,ia,ju,ao,jao,iao,ipos,job,lao,ier) 
      implicit none 
      XMPH_FLOAT a(*),ao(*) 
      MPH_INT nrow,ju(*),ia(*),ja(*),jao(*),iao(*),ipos,job,lao,ier 
c----------------------------------------------------------------------
c copies partial rows of the matrix a,ja, ia, into ao, jao, iao. 
c for each row i copies entries pointed to by ia(i) to ju(i)-1
c useful for extracting+copying lower/upper triangular parts of 
c matrices in CSR format.  
c----------------------------------------------------------------------
c on entry:
c---------
c nrow	= row dimension of the matrix 
c a,
c ja,
c ia    = input matrix in compressed sparse row format. 
c ipos  = integer. indicates the position in the array ao, jao
c         where the first element should be copied. Thus 
c         iao(1) = ipos on return. 
c job   = job indicator. if (job .ne. 1) the values are not copied 
c         (i.e., pattern only is copied in the form of arrays ja, ia).
c
c lao   = rightmost bound for arrays ao, jao. If this bound is exceeded
c         while copying, the code will exit with ier =1. 
c
c on return:
c----------
c ao,
c jao,
c iao   = output matrix containing the same data as a, ja, ia.
c ier   = 0 for normal return,
c       = 1 if size of arrays ao. jao is insufficient. 
c-----------------------------------------------------------------------
c     local variables
      integer ko, i, k 
c=======================================================================
      ier = 0 
      ko = ipos 
      iao(1) = ko
c     
      do i=1, nrow 
         do k=ia(i),ju(i)-1
            if (job .eq. 1) ao(ko) = a(k)   
            jao(ko) = ja(k)
            ko = ko+1
            if (ko .gt. lao) goto 10
         enddo
         iao(i+1) = ko 
      enddo
      return
 10    ier = 1
       return
c--------end-of-pcopmat ------------------------------------------------ 
c-----------------------------------------------------------------------
      end
c-----------------------------------------------------------------------






