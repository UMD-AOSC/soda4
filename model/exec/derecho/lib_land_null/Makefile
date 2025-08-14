# Makefile created by mkmf $Id: mkmf,v 1.1.2.1 2013/12/18 17:47:54 Niki.Zadeh Exp $ 

SRCROOT = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/

CPPDEFS = -Duse_netCDF -Duse_libMPI -DUSE_OCEAN_BGC -DENABLE_ODA -DSPMD -DLAND_BND_TRACERS -DENABLE_GDS -Duse_netCDF4

OTHERFLAGS = -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_FMS

include /glade/cheyenne/scratch/chepurin/SODA4_TEST/bin/mkmf.template.derecho


.DEFAULT:
	-echo $@ does not exist.
all: lib_land_null.a
land_model.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/land_null/land_model.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/land_null/land_model.F90
./land_model.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/land_null/land_model.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/land_null/land_model.F90 .
SRC = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/land_null/land_model.F90
OBJ = land_model.o
OFF = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/land_null/land_model.F90
clean: neat
	-rm -f .lib_land_null.a.cppdefs $(OBJ) lib_land_null.a
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
lib_land_null.a: $(OBJ)
	$(AR) $(ARFLAGS) lib_land_null.a $(OBJ)
