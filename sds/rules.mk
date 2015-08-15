# maphys dense matrix module/library
libxmph_sds_a_SOURCES = $(addprefix $(topsrcdir)/sds/, \
	xmph_sds_mod.F90 mph_sds_enum.F90 \
	xmph_sds_mumps_mod.F90 xmph_sds_mumps_type.F90 \
	xmph_sds_pastix_mod.F90 mph_sds_pastixc.c      )

libcmph_sds_a_SOURCES = $(subst xmph_,cmph_,$(libxmph_sds_a_SOURCES))
libdmph_sds_a_SOURCES = $(subst xmph_,dmph_,$(libxmph_sds_a_SOURCES))
libsmph_sds_a_SOURCES = $(subst xmph_,smph_,$(libxmph_sds_a_SOURCES))
libzmph_sds_a_SOURCES = $(subst xmph_,zmph_,$(libxmph_sds_a_SOURCES))

libcmph_sds_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libcmph_sds_a_SOURCES)))
libdmph_sds_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libdmph_sds_a_SOURCES)))
libsmph_sds_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libsmph_sds_a_SOURCES)))
libzmph_sds_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libzmph_sds_a_SOURCES)))

libcmph_sds_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libcmph_sds_a_SOURCES)))
libdmph_sds_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libdmph_sds_a_SOURCES)))
libsmph_sds_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libsmph_sds_a_SOURCES)))
libzmph_sds_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libzmph_sds_a_SOURCES)))

$(libcmph_sds_a_OBJECTS) : EXTRA_FCFLAGS=$(MUMPS_FCFLAGS) $(PASTIX_FCFLAGS)
$(libdmph_sds_a_OBJECTS) : EXTRA_FCFLAGS=$(MUMPS_FCFLAGS) $(PASTIX_FCFLAGS)
$(libsmph_sds_a_OBJECTS) : EXTRA_FCFLAGS=$(MUMPS_FCFLAGS) $(PASTIX_FCFLAGS)
$(libzmph_sds_a_OBJECTS) : EXTRA_FCFLAGS=$(MUMPS_FCFLAGS) $(PASTIX_FCFLAGS)

### 

ifneq (,$(filter c,$(ARITHS)))
include $(libcmph_sds_a_OBJECTS:.o=.d)
endif
#
ifneq (,$(filter d,$(ARITHS)))
include $(libdmph_sds_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter s,$(ARITHS)))
include $(libsmph_sds_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter z,$(ARITHS)))
include $(libzmph_sds_a_OBJECTS:.o=.d)
endif

