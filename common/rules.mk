# maphys common library (common means common to all arithmetics)
libmph_common_a_SOURCES = $(addprefix $(topsrcdir)/common/,\
	mph_common.F90 \
	mph_dbg_mod.F90 mph_env_mod.F90 \
	mph_error_mod.F90 mph_log_mod.F90 mph_mem_mod.F90 \
	mph_time_mod.F90)
libmph_common_a_OBJECTS = $(libmph_common_a_SOURCES:%.F90=%$(OBJEXT))

### 
include $(libmph_common_a_OBJECTS:.o=.d)