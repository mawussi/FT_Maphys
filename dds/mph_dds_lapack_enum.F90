!> enumerations for [CDSZ]MPH_dds_lapack_mod modules
Module MPH_dds_lapack_enum
  
  !> sizes
  Integer, parameter :: MAPHYS_LAPACK_STATE_SIZE = 6
  Integer, parameter :: MAPHYS_LAPACK_OPTION_SIZE = 2
  Integer, parameter :: MAPHYS_LAPACK_RINFO_SIZE = 2
  
  !> States
  Integer, Parameter :: LAPACK_STATE_isValid      = 1
  Integer, Parameter :: LAPACK_STATE_Factorized   = 2
  
  !> statistics
  Integer, Parameter :: LAPACK_OPTION_CompStats = 2 
  Integer, Parameter :: LAPACK_RINFO_Norm1 = 1 
  Integer, Parameter :: LAPACK_RINFO_RCond1 = 2 
  
End Module MPH_dds_lapack_enum
