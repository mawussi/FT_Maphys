c Warning: XMPH_GENFILE_COMMENT
c----------------------------------------------------------------------
c
c Partial ILUT - with Azzam modifications 
c
c----------------------------------------------------------------------
c
c This is a cleanup version of the file "xmp_pilut_commented.F".
c Please use "xmp_pilut_commented.F" if you need  
c to understand the differences between it and the original version,
c or have statistics about its execution.
c
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c Preprocessing
c----------------------------------------------------------------------

C define MPH_INT
#include "mph_defs_f.h"
#include "mph_macros_f.h"

C CONVERT2FLOAT
#if XMPH_HAVE_ARITH_C
#define CONVERT2FLOAT(x) CMPLX(x,KIND=4)       
#elif XMPH_HAVE_ARITH_D
#define CONVERT2FLOAT(x) REAL(x,KIND=8)       
#elif XMPH_HAVE_ARITH_S
#define CONVERT2FLOAT(x) REAL(x,KIND=4)       
#elif XMPH_HAVE_ARITH_Z
#define CONVERT2FLOAT(x) CMPLX(x,KIND=8)       
#endif

c----------------------------------------------------------------------- 
c PARTIAL ILUT - adapted from ILUT. Mon Aug 24 11:48:10 CDT 1998 YS. 
c-----------------------------------------------------------------------

      SUBROUTINE XMPH_pilut
     *     (n,a,ja,ia,lfil,droptol,last,lalu,alu,jlu,ju,
     *     w,jw,ierr,lfilschur,dropschur,tpilu)

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
C The initial routine is by Y. Saad for Arms. We change it to compute
C an sparse approximation of the Schur complement based on an ILU 
C factorization of the (1,1) block B.
C The ILU factorization is based on a right-looking approach. However,
C because the matrix is given in a CSR format, it should be called
C a "up-looking" variant.
C As a consequence, because we are only interested in the approximated
C Schur complement, the entries computed in the L part are not kept.
C In addition, we used two different thresholds, the first "droptol" is 
C the regular ILU threshold in the U and U1 matrices. The second one is 
C "dropschur", that is used for the droppting in the approximated Schur.
C Similarly and extra fil parameter has been added to enable more fill
C in the Schur complement than in the U and U1 parts.
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

         if((j2-j1+1).eq.0) Goto 993

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

            lenu = len+1
            len = min0(lenu,lfil)
         else ! Schur part U
            len = 0
            do k=1, lenu-1
               if (ABS(w(ibound+k)) .gt. dropschur*tnorm) then 
                  len = len+1
                  w(ibound+len) = w(ibound+k) 
                  jw(ibound+len) = jw(ibound+k) 
               endif
            enddo
            lenu = len+1
            len = min0(lenu,lfilschur)
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






