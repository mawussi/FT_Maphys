# maphys dense matrix module/library
xmph_test_all_SOURCES = $(addprefix $(topsrcdir)/unittest/, \
	xmph_test_all.F90 \
	xmph_test_fixture_mod.F90 \
	xmph_test_dense_matrix_mod.F90 \
	xmph_test_sparse_matrix_mod.F90 \
	xmph_test_dds_mod.F90 \
	xmph_test_sds_mod.F90 \
	xmph_test_ilu_mod.F90 \
	xmph_test_part_mod.F90 \
	test_utils_mod.F90)

xmph_test_all_LDADD = $(addprefix $(buildlibdir)/,\
	libmaphys.a libxmph_toolkit.a libslatec.a )

cmph_test_all_SOURCES = $(subst xmph_,cmph_,$(xmph_test_all_SOURCES))
dmph_test_all_SOURCES = $(subst xmph_,dmph_,$(xmph_test_all_SOURCES))
smph_test_all_SOURCES = $(subst xmph_,smph_,$(xmph_test_all_SOURCES))
zmph_test_all_SOURCES = $(subst xmph_,zmph_,$(xmph_test_all_SOURCES))

cmph_test_all_LDADD = $(subst xmph_,cmph_,$(xmph_test_all_LDADD))
dmph_test_all_LDADD = $(subst xmph_,dmph_,$(xmph_test_all_LDADD))
smph_test_all_LDADD = $(subst xmph_,smph_,$(xmph_test_all_LDADD))
zmph_test_all_LDADD = $(subst xmph_,zmph_,$(xmph_test_all_LDADD))

cmph_test_all_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(cmph_test_all_SOURCES)))
dmph_test_all_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(dmph_test_all_SOURCES)))
smph_test_all_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(smph_test_all_SOURCES)))
zmph_test_all_OBJECTS += $(patsubst %.F90,%.o,$(filter %.F90,$(zmph_test_all_SOURCES)))

cmph_test_all_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(cmph_test_all_SOURCES)))
dmph_test_all_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(dmph_test_all_SOURCES)))
smph_test_all_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(smph_test_all_SOURCES)))
zmph_test_all_OBJECTS += $(patsubst %.c,%.o,$(filter %.c,$(zmph_test_all_SOURCES)))

$(topsrcdir)/unittest/cmph_test_part_mod.o : EXTRA_FCFLAGS= $(METIS_FCFLAGS) $(SCOTCH_FCFLAGS) 
$(topsrcdir)/unittest/dmph_test_part_mod.o : EXTRA_FCFLAGS= $(METIS_FCFLAGS) $(SCOTCH_FCFLAGS) 
$(topsrcdir)/unittest/smph_test_part_mod.o : EXTRA_FCFLAGS= $(METIS_FCFLAGS) $(SCOTCH_FCFLAGS) 
$(topsrcdir)/unittest/zmph_test_part_mod.o : EXTRA_FCFLAGS= $(METIS_FCFLAGS) $(SCOTCH_FCFLAGS) 

###

ifneq (,$(filter c,$(ARITHS)))
include $(cmph_test_all_OBJECTS:.o=.d)
endif
#
ifneq (,$(filter d,$(ARITHS)))
include $(dmph_test_all_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter s,$(ARITHS)))
include $(smph_test_all_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter z,$(ARITHS)))
include $(zmph_test_all_OBJECTS:.o=.d)
endif

