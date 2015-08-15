# maphys dense matrix module/library
libxmph_dds_a_SOURCES = $(addprefix $(topsrcdir)/dds/, \
	xmph_dds_mod.F90 \
	mph_dds_enum.F90 \
	xmph_dds_lapack_mod.F90 \
	mph_dds_lapack_enum.F90 )

libcmph_dds_a_SOURCES = $(subst xmph_,cmph_,$(libxmph_dds_a_SOURCES))
libdmph_dds_a_SOURCES = $(subst xmph_,dmph_,$(libxmph_dds_a_SOURCES))
libsmph_dds_a_SOURCES = $(subst xmph_,smph_,$(libxmph_dds_a_SOURCES))
libzmph_dds_a_SOURCES = $(subst xmph_,zmph_,$(libxmph_dds_a_SOURCES))

libcmph_dds_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libcmph_dds_a_SOURCES)))
libdmph_dds_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libdmph_dds_a_SOURCES)))
libsmph_dds_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libsmph_dds_a_SOURCES)))
libzmph_dds_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libzmph_dds_a_SOURCES)))

libcmph_dds_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libcmph_dds_a_SOURCES)))
libdmph_dds_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libdmph_dds_a_SOURCES)))
libsmph_dds_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libsmph_dds_a_SOURCES)))
libzmph_dds_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libzmph_dds_a_SOURCES)))

$(libcmph_dds_a_OBJECTS) : EXTRA_FCFLAGS=$(DALGEBRA_FCFLAGS)
$(libdmph_dds_a_OBJECTS) : EXTRA_FCFLAGS=$(DALGEBRA_FCFLAGS)
$(libsmph_dds_a_OBJECTS) : EXTRA_FCFLAGS=$(DALGEBRA_FCFLAGS)
$(libzmph_dds_a_OBJECTS) : EXTRA_FCFLAGS=$(DALGEBRA_FCFLAGS)

### 

ifneq (,$(filter c,$(ARITHS)))
include $(libcmph_dds_a_OBJECTS:.o=.d)
endif
#
ifneq (,$(filter d,$(ARITHS)))
include $(libdmph_dds_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter s,$(ARITHS)))
include $(libsmph_dds_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter z,$(ARITHS)))
include $(libzmph_dds_a_OBJECTS:.o=.d)
endif

