#
#  set compilation flags for MaPHyS
#
#  FCFLAGS   -> compil F90 code
#  CFLAGS    -> compil C   code 
#  FDFLAGS   -> link objects

#-------------------------------------------------------------------------------
# [2.0] PARTITIONNERS
#-------------------------------------------------------------------------------

#---------------------------------------------------------------------
# [2.1] METIS - http://glaros.dtc.umn.edu/gkhome/views/metis
#---------------------------------------------------------------------

ifdef METIS_DIR
METIS_def = -DHAVE_METIS
METIS_inc = -I$(METIS_DIR)/Lib
METIS_lib = -L$(METIS_DIR)/ -lmetis
endif

#---------------------------------------------------------------------
# [2.2] SCOTCH - http://gforge.inria.fr/projects/scotch/
#---------------------------------------------------------------------

ifdef SCOTCH_DIR
SCOTCH_def ?= -DHAVE_SCOTCH
SCOTCH_inc ?= -I$(SCOTCH_DIR)/include
SCOTCH_lib ?= -L$(SCOTCH_DIR)/lib -lscotch -lscotcherrexit 
endif

#---------------------------------------------------------------------
# [2.3] SUMMARY
#---------------------------------------------------------------------

PART_def = $(METIS_def) $(SCOTCH_def)
PART_inc = $(METIS_inc) $(SCOTCH_inc)
PART_lib = $(METIS_lib) $(SCOTCH_lib) $(AZZLIBPAR)

#-------------------------------------------------------------------------------
# [3.0] DIRECT_SOLVERS
#-------------------------------------------------------------------------------

#---------------------------------------------------------------------
# [3.1] MUMPS  - http://mumps.enseeiht.fr/
#       PaSTiX - http://gforge.inria.fr/projects/pastix/
#---------------------------------------------------------------------

include $(topdir)/config/mumps.mk
include $(topdir)/config/pastix.mk

#---------------------------------------------------------------------
# [3.2] HWLOC - http://www.open-mpi.org/projects/hwloc/
#---------------------------------------------------------------------

include $(topdir)/config/hwloc.mk

#---------------------------------------------------------------------
# [3.3] Sum up
#---------------------------------------------------------------------

DSLV_def = $(MUMPS_def) $(PASTIX_def) $(HWLOC_def)
DSLV_inc = $(MUMPS_inc) $(PASTIX_inc) $(HWLOC_inc)
DSLV_lib = $(MUMPS_lib) $(PASTIX_lib) $(HWLOC_lib)

#-------------------------------------------------------------------------------
# [4.0] MPI + Mathematical libraries (scalapack, blacs, lapack, blas)
#-------------------------------------------------------------------------------

MPI_inc    = 
MPI_lib    = 

MAPHYS_MATH_lib = $(SCALAP) $(BLACS) $(BLAS) $(LAPACK) $(BLAS) $(LIBMPI) $(MATHLIB)

#-------------------------------------------------------------------------------
# [5.0] SUMMARY of all headers and libraries
#-------------------------------------------------------------------------------

MAPHYS_def := $(DSLV_def) $(PART_def)

MAPHYS_inc := -I$(topdir)/include -I$(moddir)
MAPHYS_inc += $(DSLV_inc) $(PART_inc) $(MPI_inc) 

#MAPHYS_lib := $(MPHS_lib) $(DSLV_lib) $(PART_lib) $(MAPHYS_MATH_lib)
MAPHYS_lib := $(DSLV_lib) $(PART_lib) $(MAPHYS_MATH_lib)

CFLAGS  += $(MAPHYS_def) $(MAPHYS_inc)
FCFLAGS += $(MAPHYS_def) $(MAPHYS_inc)
FDFLAGS += $(MAPHYS_lib)
