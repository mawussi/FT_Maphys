#---- PASTIX ----
ifdef PASTIX_DIR
PASTIX_def   = -DHAVE_LIBPASTIX
PASTIX_inc   = -I$(PASTIX_DIR)
PASTIX_lib  := -L$(PASTIX_DIR) -lpastix -lrt
endif
