# @file Makefile
#
# MaPHyS is a software package provided by INRIA.
#
# Build libmaphys.a or lib[cdsz]maphys.a, and copy its dependencies
#
# @author Yohan Lee-tin-yien
#
###

libcmaphys_a_LIBADD = $(libmph_common_a_OBJECTS) \
	$(libcmph_dm_a_OBJECTS) $(libcmph_sm_a_OBJECTS) \
	$(libcmph_dds_a_OBJECTS) $(libcmph_sds_a_OBJECTS) \
	$(libcmph_ilu_a_OBJECTS) \
	$(libcmph_part_a_OBJECTS) $(libcmph_hybrid_a_OBJECTS) 

libdmaphys_a_LIBADD = $(libmph_common_a_OBJECTS) \
	$(libdmph_dm_a_OBJECTS) $(libdmph_sm_a_OBJECTS) \
	$(libdmph_dds_a_OBJECTS) $(libdmph_sds_a_OBJECTS) \
	$(libdmph_ilu_a_OBJECTS) \
	$(libdmph_part_a_OBJECTS) $(libdmph_hybrid_a_OBJECTS) 

libsmaphys_a_LIBADD = $(libmph_common_a_OBJECTS) \
	$(libsmph_dm_a_OBJECTS) $(libsmph_sm_a_OBJECTS) \
	$(libsmph_dds_a_OBJECTS) $(libsmph_sds_a_OBJECTS) \
	$(libsmph_ilu_a_OBJECTS) \
	$(libsmph_part_a_OBJECTS) $(libsmph_hybrid_a_OBJECTS) 

libzmaphys_a_LIBADD = $(libmph_common_a_OBJECTS) \
	$(libzmph_dm_a_OBJECTS) $(libzmph_sm_a_OBJECTS) \
	$(libzmph_dds_a_OBJECTS) $(libzmph_sds_a_OBJECTS) \
	$(libzmph_ilu_a_OBJECTS) \
	$(libzmph_part_a_OBJECTS) $(libzmph_hybrid_a_OBJECTS) 
##
libmaphys_a_LIBADD = 

ifneq (,$(filter c,$(ARITHS)))
libmaphys_a_LIBADD += $(libcmaphys_a_LIBADD)
endif

ifneq (,$(filter d,$(ARITHS)))
libmaphys_a_LIBADD += $(libdmaphys_a_LIBADD)
endif

ifneq (,$(filter z,$(ARITHS)))
libmaphys_a_LIBADD += $(libsmaphys_a_LIBADD)
endif

ifneq (,$(filter z,$(ARITHS)))
libmaphys_a_LIBADD += $(libzmaphys_a_LIBADD)
endif
