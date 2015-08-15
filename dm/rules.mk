# maphys dense matrix module/library

libxmph_dm_a_SOURCES = $(topsrcdir)/dm/xmph_dense_matrix_mod.F90 $(topsrcdir)/dm/xmph_dense_matrix_type.F90

libcmph_dm_a_SOURCES = $(subst xmph_,cmph_,$(libxmph_dm_a_SOURCES))
libdmph_dm_a_SOURCES = $(subst xmph_,dmph_,$(libxmph_dm_a_SOURCES))
libsmph_dm_a_SOURCES = $(subst xmph_,smph_,$(libxmph_dm_a_SOURCES))
libzmph_dm_a_SOURCES = $(subst xmph_,zmph_,$(libxmph_dm_a_SOURCES))

libmph_dm_a_SOURCES = $(topsrcdir)/dm/mph_dense_matrix_mod.F90                 \
	$(sort $(libcmph_dm_a_SOURCES) $(libdmph_dm_a_SOURCES) \
	$(libsmph_dm_a_SOURCES) $(libzmph_dm_a_SOURCES))

libcmph_dm_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libcmph_dm_a_SOURCES)))
libdmph_dm_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libdmph_dm_a_SOURCES)))
libsmph_dm_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libsmph_dm_a_SOURCES)))
libzmph_dm_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libzmph_dm_a_SOURCES)))

libmph_dm_a_OBJECTS  += $(patsubst %.F90,%.o,$(filter %.F90,$(libmph_dm_a_SOURCES))) 

### 

ifneq (,$(filter c,$(ARITHS)))
include $(libcmph_dm_a_OBJECTS:.o=.d)
endif
#
ifneq (,$(filter d,$(ARITHS)))
include $(libdmph_dm_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter s,$(ARITHS)))
include $(libsmph_dm_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter z,$(ARITHS)))
include $(libzmph_dm_a_OBJECTS:.o=.d)
endif


