! Warning: XMPH_GENFILE_COMMENT
! Module including the mumps type definition.
!
! This is a trick to allow structure renaming.
Module XMPH_sds_mumps_type

  !* No implicit typing *!
  Implicit None

  !* Type definition *! 

  ! define XMPH_ARITHmumps_struc

#if HAVE_LIBMUMPS

  Include 'XMPH_ARITHmumps_struc.h'

#else

  Type XMPH_ARITHmumps_struc; Sequence
     Integer :: dummy
  End type XMPH_ARITHmumps_struc

#endif

End Module XMPH_sds_mumps_type
