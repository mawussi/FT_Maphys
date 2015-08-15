! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"

!> module containing the routines to partition/distribute the system. 
!!
!! @author Azzam Haidar
!! @author Yohan Lee-tin-yien
!!
Module XMPH_part_mod
  
  !* Imported types *!
  Use MPH_part_type

  !* Imported routines *!
  Use MPH_domain_mod, Only : MPH_domain_write

  Use XMPH_part_ordermatrix_mod, Only : XMPH_PART_OrderGlobalMatrix
  Use XMPH_part_builddomains_mod, Only : XMPH_PART_Build_Domains
  Use XMPH_part_distmatrix_mod, Only : &
       XMPH_PART_DistGlobalMatrix, XMPH_PART_Compute_interface_weight
  Use XMPH_part_distrhs_mod, Only : XMPH_PART_DistRHS
  Use XMPH_part_collectsol_mod, Only : XMPH_PART_CollectSol


End Module XMPH_part_mod
