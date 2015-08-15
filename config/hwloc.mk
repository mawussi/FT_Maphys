#---- HWLOC ----

ifeq (-DHAVE_HWLOC,$(findstring -DHAVE_HWLOC,$(FCFLAGS)))

ifdef HWLOC_prefix
HWLOC_inc ?= -I$(HWLOC_prefix)/include
HWLOC_lib ?= -L$(HWLOC_prefix)/lib -lhwloc
else
HWLOC_inc ?= 
HWLOC_lib ?= -lhwloc
endif


endif