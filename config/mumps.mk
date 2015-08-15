#---- MUMPS ----
ifdef MUMPS_DIR

MUMPS_def       += -DHAVE_LIBMUMPS
MUMPS_inc 	= -I$(MUMPS_DIR)/include
MUMPS_lib   	= -L$(MUMPS_DIR)/lib \
	          -l$(ARITH)mumps -lmumps_common \
		  -L$(MUMPS_DIR)/PORD/lib -lpord
endif