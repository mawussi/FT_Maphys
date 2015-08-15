#include "pastix_fortran.h"

#if   defined(MAPHYS_ARITH_d)
#define x_pastix_fortran                     d_pastix_fortran
#define x_pastix_fortran_checkmatrix         d_pastix_fortran_checkmatrix
#define x_pastix_fortran_checkmatrix_end     d_pastix_fortran_checkmatrix_end
#define x_pastix_fortran_setschurunknownlist d_pastix_fortran_setschurunknownlist 
#define x_pastix_fortran_setschurarray       d_pastix_fortran_setschurarray
#define x_pastix_fortran_getSchur            d_pastix_fortran_getSchur            
#define x_pastix_fortran_setbindtab          d_pastix_fortran_setbindtab 
#define x_pastix_fortran_getschurlocalnodenbr d_pastix_fortran_getschurlocalnodenbr

#elif defined(MAPHYS_ARITH_s)
#define x_pastix_fortran                     s_pastix_fortran
#define x_pastix_fortran_checkmatrix         s_pastix_fortran_checkmatrix
#define x_pastix_fortran_checkmatrix_end     s_pastix_fortran_checkmatrix_end
#define x_pastix_fortran_setschurunknownlist s_pastix_fortran_setschurunknownlist 
#define x_pastix_fortran_setschurarray       s_pastix_fortran_setschurarray
#define x_pastix_fortran_getSchur            s_pastix_fortran_getSchur            
#define x_pastix_fortran_setbindtab          s_pastix_fortran_setbindtab 
#define x_pastix_fortran_getschurlocalnodenbr s_pastix_fortran_getschurlocalnodenbr

#elif defined(MAPHYS_ARITH_z)
#define x_pastix_fortran                     z_pastix_fortran
#define x_pastix_fortran_checkmatrix         z_pastix_fortran_checkmatrix
#define x_pastix_fortran_checkmatrix_end     z_pastix_fortran_checkmatrix_end
#define x_pastix_fortran_setschurunknownlist z_pastix_fortran_setschurunknownlist 
#define x_pastix_fortran_setschurarray       z_pastix_fortran_setschurarray
#define x_pastix_fortran_getSchur            z_pastix_fortran_getSchur   
#define x_pastix_fortran_setbindtab          z_pastix_fortran_setbindtab          
#define x_pastix_fortran_getschurlocalnodenbr z_pastix_fortran_getschurlocalnodenbr

#elif defined(MAPHYS_ARITH_c)

#define x_pastix_fortran                     c_pastix_fortran
#define x_pastix_fortran_checkmatrix         c_pastix_fortran_checkmatrix
#define x_pastix_fortran_checkmatrix_end     c_pastix_fortran_checkmatrix_end
#define x_pastix_fortran_setschurunknownlist c_pastix_fortran_setschurunknownlist 
#define x_pastix_fortran_setschurarray       c_pastix_fortran_setschurarray
#define x_pastix_fortran_getSchur            c_pastix_fortran_getSchur   
#define x_pastix_fortran_setbindtab          c_pastix_fortran_setbindtab          
#define x_pastix_fortran_getschurlocalnodenbr c_pastix_fortran_getschurlocalnodenbr

#else
#error "wrong MAPHYS_ARITH value"
#endif

! forgotten index in fortran api
#ifndef IPARM_SCHUR
#define IPARM_SCHUR 28
#endif

