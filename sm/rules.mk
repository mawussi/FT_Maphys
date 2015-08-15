# maphys dense matrix module/library
libxmph_sm_a_SOURCES = $(addprefix $(topsrcdir)/sm/, \
	xmph_sparse_matrix_mod.F90 \
	xmph_sparse_matrix_type.F90 \
	mph_qsortc.c \
	mph_sparse_matrix_enum.F90  )

libcmph_sm_a_SOURCES = $(subst xmph_,cmph_,$(libxmph_sm_a_SOURCES))
libdmph_sm_a_SOURCES = $(subst xmph_,dmph_,$(libxmph_sm_a_SOURCES))
libsmph_sm_a_SOURCES = $(subst xmph_,smph_,$(libxmph_sm_a_SOURCES))
libzmph_sm_a_SOURCES = $(subst xmph_,zmph_,$(libxmph_sm_a_SOURCES))

libcmph_sm_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libcmph_sm_a_SOURCES)))
libdmph_sm_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libdmph_sm_a_SOURCES)))
libsmph_sm_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libsmph_sm_a_SOURCES)))
libzmph_sm_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libzmph_sm_a_SOURCES)))

libcmph_sm_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libcmph_sm_a_SOURCES)))
libdmph_sm_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libdmph_sm_a_SOURCES)))
libsmph_sm_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libsmph_sm_a_SOURCES)))
libzmph_sm_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libzmph_sm_a_SOURCES)))

### 

ifneq (,$(filter c,$(ARITHS)))
include $(libcmph_sm_a_OBJECTS:.o=.d)
endif
#
ifneq (,$(filter d,$(ARITHS)))
include $(libdmph_sm_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter s,$(ARITHS)))
include $(libsmph_sm_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter z,$(ARITHS)))
include $(libzmph_sm_a_OBJECTS:.o=.d)
endif

