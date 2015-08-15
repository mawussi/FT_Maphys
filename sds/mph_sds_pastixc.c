/* C functions for the wrapper sds_pastix */

#if PASTIX_HAS_MEMORY_USAGE

unsigned long memAllocGetCurrent ();
void mph_sds_pastix_getmemusage(long *memusage )
{  *memusage = (long) memAllocGetCurrent (); }

#else

void mph_sds_pastix_getmemusage(long *memusage )
{  *memusage = -1L ;}

#endif


/* Fortran name mangling */

void MPH_SDS_PASTIX_GETMEMUSAGE(long *memusage )
{
  mph_sds_pastix_getmemusage(memusage);
}
void mph_sds_pastix_getmemusage_(long *memusage )
{
  mph_sds_pastix_getmemusage(memusage);
}
void mph_sds_pastix_getmemusage__(long *memusage )
{
  mph_sds_pastix_getmemusage(memusage);
}

