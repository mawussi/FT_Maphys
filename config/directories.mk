### BUILDING DIRECTORIES ###
# directories for libraries
libdir  := $(libdir)

# directories for maphys
# object files
dirs    := $(foreach o,$(OBJS)  , $(dir $(o)) )
blddirs := $(sort $(dirs))

# building rules
$(blddirs): ;	mkdir -p $@
$(libdir) : ; 	mkdir -p $@
$(moddir) : ; 	mkdir -p $@
