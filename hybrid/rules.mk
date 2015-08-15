# maphys dense matrix module/library

libxmph_hybrid_a_SOURCES = $(addprefix $(topsrcdir)/hybrid/, \
	xmph_dls_mod.F90 xmph_sls_mod.F90 \
	mph_maphys_enum.F90 xmph_maphys_type.F90 \
	xmph_state_mod.F90 xmph_state_update_mod.F90 xmph_state_print_mod.F90 \
	mph_strat_mod.F90 mph_strat_utils_mod.F90 \
	xmph_distsystem_mod.F90 \
	xmph_domainsls_mod.F90 xmph_schur_assemble_mod.F90 \
	xmph_pilut_mod.F90 xmph_pilut.F \
	xmph_schur_aux_mod.F90 xmph_schur_mod.F90 \
	xmph_pcd_mod.F90 xmph_maphys_aux_mod.F90 xmph_maphys_mod.F90 xmph_fault_mod.F90 )

libcmph_hybrid_a_SOURCES = $(subst xmph_,cmph_,$(libxmph_hybrid_a_SOURCES))
libdmph_hybrid_a_SOURCES = $(subst xmph_,dmph_,$(libxmph_hybrid_a_SOURCES))
libsmph_hybrid_a_SOURCES = $(subst xmph_,smph_,$(libxmph_hybrid_a_SOURCES))
libzmph_hybrid_a_SOURCES = $(subst xmph_,zmph_,$(libxmph_hybrid_a_SOURCES))

libcmph_hybrid_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libcmph_hybrid_a_SOURCES)))
libdmph_hybrid_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libdmph_hybrid_a_SOURCES)))
libsmph_hybrid_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libsmph_hybrid_a_SOURCES)))
libzmph_hybrid_a_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(libzmph_hybrid_a_SOURCES)))

libcmph_hybrid_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libcmph_hybrid_a_SOURCES)))
libdmph_hybrid_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libdmph_hybrid_a_SOURCES)))
libsmph_hybrid_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libsmph_hybrid_a_SOURCES)))
libzmph_hybrid_a_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(libzmph_hybrid_a_SOURCES)))

libcmph_hybrid_a_OBJECTS += $(patsubst %.F,%.o,$(filter %.F,$(libcmph_hybrid_a_SOURCES)))
libdmph_hybrid_a_OBJECTS += $(patsubst %.F,%.o,$(filter %.F,$(libdmph_hybrid_a_SOURCES)))
libsmph_hybrid_a_OBJECTS += $(patsubst %.F,%.o,$(filter %.F,$(libsmph_hybrid_a_SOURCES)))
libzmph_hybrid_a_OBJECTS += $(patsubst %.F,%.o,$(filter %.F,$(libzmph_hybrid_a_SOURCES)))

### 

ifneq (,$(filter c,$(ARITHS)))
include $(libcmph_hybrid_a_OBJECTS:.o=.d)
endif
#
ifneq (,$(filter d,$(ARITHS)))
include $(libdmph_hybrid_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter s,$(ARITHS)))
include $(libsmph_hybrid_a_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter z,$(ARITHS)))
include $(libzmph_hybrid_a_OBJECTS:.o=.d)
endif

