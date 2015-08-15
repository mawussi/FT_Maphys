# maphys dense matrix module/library
libxmph_part_a_SOURCES = $(addprefix $(topsrcdir)/part/, \
	mph_part_type.F90 \
	mph_part_scotch_mod.F90 \
	mph_domain_mod.F90 \
	mph_rhs_partition_mod.F90 \
	xmph_part_mod.F90 \
	xmph_part_ordermatrix_mod.F90 \
	xmph_part_builddomains_mod.F90 \
	xmph_part_distmatrix_mod.F90 \
	xmph_part_distrhs_mod.F90 \
	xmph_part_collectsol_mod.F90 )

# Special case for METIS.
# Note: this is uneccessary in for mph_part_scotch_mod.F90 which already handles internally -DHAVE_LIBSCOTCH
ifneq ($(filter $(METIS_FCFLAGS),-DHAVE_LIBMETIS),)
libxmph_part_a_SOURCES += $(addprefix $(topsrcdir)/part/, \
	mph_part_metis.c mph_part_metis_fortran.c \
	mph_part_metis_proto.h mph_part_metis_struct.h )
endif


libcmph_part_a_SOURCES = $(subst xmph_,cmph_,$(libxmph_part_a_SOURCES))
libdmph_part_a_SOURCES = $(subst xmph_,dmph_,$(libxmph_part_a_SOURCES))
libsmph_part_a_SOURCES = $(subst xmph_,smph_,$(libxmph_part_a_SOURCES))
libzmph_part_a_SOURCES = $(subst xmph_,zmph_,$(libxmph_part_a_SOURCES))

libcmph_part_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libcmph_part_a_SOURCES)))
libdmph_part_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libdmph_part_a_SOURCES)))
libsmph_part_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libsmph_part_a_SOURCES)))
libzmph_part_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libzmph_part_a_SOURCES)))

libcmph_part_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libcmph_part_a_SOURCES)))
libdmph_part_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libdmph_part_a_SOURCES)))
libsmph_part_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libsmph_part_a_SOURCES)))
libzmph_part_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libzmph_part_a_SOURCES)))

# Warning : order in EXTRA_*CFLAGS defined below are important.
#
# Since we call specific subfunctions of libmetis.a.
# we must include metis.h from METIS before the one defined by SCOTCH.
#
#
$(libcmph_part_a_OBJECTS) : EXTRA_FCFLAGS= $(METIS_FCFLAGS) $(SCOTCH_FCFLAGS) 
$(libdmph_part_a_OBJECTS) : EXTRA_FCFLAGS= $(METIS_FCFLAGS) $(SCOTCH_FCFLAGS) 
$(libsmph_part_a_OBJECTS) : EXTRA_FCFLAGS= $(METIS_FCFLAGS) $(SCOTCH_FCFLAGS) 
$(libzmph_part_a_OBJECTS) : EXTRA_FCFLAGS= $(METIS_FCFLAGS) $(SCOTCH_FCFLAGS) 

$(libcmph_part_a_OBJECTS) : EXTRA_CFLAGS= $(METIS_CFLAGS) $(SCOTCH_CFLAGS)
$(libdmph_part_a_OBJECTS) : EXTRA_CFLAGS= $(METIS_CFLAGS) $(SCOTCH_CFLAGS)
$(libsmph_part_a_OBJECTS) : EXTRA_CFLAGS= $(METIS_CFLAGS) $(SCOTCH_CFLAGS)
$(libzmph_part_a_OBJECTS) : EXTRA_CFLAGS= $(METIS_CFLAGS) $(SCOTCH_CFLAGS)

### 

ifneq (,$(filter c,$(ARITHS)))
include $(libcmph_part_a_OBJECTS:.o=.d)
endif
#
ifneq (,$(filter d,$(ARITHS)))
include $(libdmph_part_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter s,$(ARITHS)))
include $(libsmph_part_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter z,$(ARITHS)))
include $(libzmph_part_a_OBJECTS:.o=.d)
endif

