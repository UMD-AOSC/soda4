# Makefile created by mkmf $Id: mkmf,v 1.1.2.1 2013/12/18 17:47:54 Niki.Zadeh Exp $ 

SRCROOT = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/

CPPDEFS = -Duse_netCDF -Duse_libMPI -DUSE_OCEAN_BGC -DENABLE_ODA -DSPMD -DLAND_BND_TRACERS -DENABLE_GDS -Duse_netCDF4


include /glade/cheyenne/scratch/chepurin/SODA4_TEST/bin/mkmf.template.derecho


.DEFAULT:
	-echo $@ does not exist.
all: lib_version.a
version.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/version/version.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/version/version.F90
./version.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/version/version.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/version/version.F90 .
SRC = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/version/version.F90
OBJ = version.o
OFF = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/version/version.F90
clean: neat
	-rm -f .lib_version.a.cppdefs $(OBJ) lib_version.a
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
lib_version.a: $(OBJ)
	$(AR) $(ARFLAGS) lib_version.a $(OBJ)
