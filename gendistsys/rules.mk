# testing partitioner

xmph_testdistsys_SOURCES := $(addprefix $(topsrcdir)/gendistsys/,\
	xmph_testdistsys.F90 xmph_gendistsys_mod.F90 \
	coordiscrete3D.f coordset_coef.f \
	wcoord_constant.f wcoord_hd.f wcoord_regions.f \
	getsizeInterf3.f setupInterf3.f )

xmph_testdistsys_LDADD = $(addprefix $(buildlibdir)/,\
	libxmph_toolkit.a libmaphys.a libslatec.a libpackcg.a libpackgmres.a libpackfgmres.a )

cmph_testdistsys_SOURCES := $(subst xmph_,cmph_,$(xmph_testdistsys_SOURCES))
dmph_testdistsys_SOURCES := $(subst xmph_,dmph_,$(xmph_testdistsys_SOURCES))
smph_testdistsys_SOURCES := $(subst xmph_,smph_,$(xmph_testdistsys_SOURCES))
zmph_testdistsys_SOURCES := $(subst xmph_,zmph_,$(xmph_testdistsys_SOURCES))

cmph_testdistsys_LDADD := $(subst xmph_,cmph_,$(xmph_testdistsys_LDADD))
dmph_testdistsys_LDADD := $(subst xmph_,dmph_,$(xmph_testdistsys_LDADD))
smph_testdistsys_LDADD := $(subst xmph_,smph_,$(xmph_testdistsys_LDADD))
zmph_testdistsys_LDADD := $(subst xmph_,zmph_,$(xmph_testdistsys_LDADD))

cmph_testdistsys_OBJECTS := $(patsubst %.F90,%.o,$(filter %.F90,$(cmph_testdistsys_SOURCES)))
dmph_testdistsys_OBJECTS := $(patsubst %.F90,%.o,$(filter %.F90,$(dmph_testdistsys_SOURCES)))
smph_testdistsys_OBJECTS := $(patsubst %.F90,%.o,$(filter %.F90,$(smph_testdistsys_SOURCES)))
zmph_testdistsys_OBJECTS := $(patsubst %.F90,%.o,$(filter %.F90,$(zmph_testdistsys_SOURCES)))

cmph_testdistsys_OBJECTS += $(patsubst %.f,%.o,$(filter %.f,$(cmph_testdistsys_SOURCES)))
dmph_testdistsys_OBJECTS += $(patsubst %.f,%.o,$(filter %.f,$(dmph_testdistsys_SOURCES)))
smph_testdistsys_OBJECTS += $(patsubst %.f,%.o,$(filter %.f,$(smph_testdistsys_SOURCES)))
zmph_testdistsys_OBJECTS += $(patsubst %.f,%.o,$(filter %.f,$(zmph_testdistsys_SOURCES)))

### 
ifneq (,$(filter c,$(ARITHS)))
include $(cmph_testdistsys_OBJECTS:.o=.d)
endif
#
ifneq (,$(filter d,$(ARITHS)))
include $(dmph_testdistsys_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter s,$(ARITHS)))
include $(smph_testdistsys_OBJECTS:.o=.d)
endif
# 
ifneq (,$(filter z,$(ARITHS)))
include $(zmph_testdistsys_OBJECTS:.o=.d)
endif

ifneq (,$(findstring YES,$(KEEP_GENERATED_SOURCES)))
.SECONDARY : $(cmph_testdistsys_SOURCES)
.SECONDARY : $(dmph_testdistsys_SOURCES)
.SECONDARY : $(smph_testdistsys_SOURCES)
.SECONDARY : $(zmph_testdistsys_SOURCES)
endif