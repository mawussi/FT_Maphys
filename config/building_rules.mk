#
# Define default building rules for maphys
#
# need for out-of-source : 
#         bd = building directory
#         sd = source   directory


#-------------------------------------------------------------------------------
# [1.0] Build directory
#-------------------------------------------------------------------------------
%/: 
	mkdir -p $@

#-------------------------------------------------------------------------------
# [1.0] Build objects from C code
#-------------------------------------------------------------------------------

# same directory rules
%.o : %.c
	$(CC) $(CFLAGS)  -I$(@D)  -c $< -o $@

# different directory rules
$(bd)/%.o : $(sd)/%.c
	$(CC) $(CFLAGS)  -c $< -o $@

#-------------------------------------------------------------------------------
# [2.0] Build objects from Fortran code
#-------------------------------------------------------------------------------

# F77
$(bd)/%.o : $(sd)/%.f
	$(FC) $(FCFLAGS) -I$(@D) -c $< -o $@	

$(bd)/%.o : $(sd)/%.F
	$(FC) $(FCFLAGS) -I$(@D) -c $< -o $@	


# F90
$(bd)/%.o : $(sd)/%.F90
	cd $(moddir) && $(FC) $(FCFLAGS) -I$(@D) -c $< -o $@	

$(bd)/%.o : $(sd)/%.f90
	cd $(moddir) && $(FC) $(FCFLAGS) -I$(@D) -c $< -o $@	

$(bd)/%.f90 : $(sd)/%.F90 
	$(CPP) $(FCFLAGS) -I$(@D) $<  $@	


