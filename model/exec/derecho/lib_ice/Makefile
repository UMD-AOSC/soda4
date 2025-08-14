# Makefile created by mkmf $Id: mkmf,v 1.1.2.1 2013/12/18 17:47:54 Niki.Zadeh Exp $ 

SRCROOT = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/

CPPDEFS = -Duse_netCDF -Duse_libMPI -DUSE_OCEAN_BGC -DENABLE_ODA -DSPMD -DLAND_BND_TRACERS -DENABLE_GDS -Duse_netCDF4

OTHERFLAGS = -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_FMS

include /glade/cheyenne/scratch/chepurin/SODA4_TEST/bin/mkmf.template.derecho


.DEFAULT:
	-echo $@ does not exist.
all: lib_ice.a
ice_albedo.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_param/ice_albedo.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_param/ice_albedo.F90
ice_bergs.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_bergs.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_bergs.F90
ice_dyn.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_dyn.F90 ice_grid.o ice_thm.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_dyn.F90
ice_grid.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_grid.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_grid.F90
ice_model.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_model.F90 ocean_albedo.o ocean_rough.o ice_type.o ice_grid.o ice_spec.o ice_thm.o ice_dyn.o ice_bergs.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_model.F90
ice_spec.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_spec.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_spec.F90
ice_thm.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_thm.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_thm.F90
ice_type.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_type.F90 ice_grid.o ice_thm.o ice_dyn.o ice_bergs.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_type.F90
mask.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/mask.F90 rot.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/mask.F90
ocean_albedo.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_param/ocean_albedo.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_param/ocean_albedo.F90
ocean_rough.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_param/ocean_rough.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_param/ocean_rough.F90
rot.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/rot.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/rot.F90
./ocean_rough.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_param/ocean_rough.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_param/ocean_rough.F90 .
./ice_albedo.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_param/ice_albedo.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_param/ice_albedo.F90 .
./ice_bergs.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_bergs.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_bergs.F90 .
./ocean_albedo.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_param/ocean_albedo.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_param/ocean_albedo.F90 .
./ice_thm.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_thm.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_thm.F90 .
./ice_spec.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_spec.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_spec.F90 .
./ice_model.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_model.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_model.F90 .
./ice_type.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_type.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_type.F90 .
./ice_grid.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_grid.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_grid.F90 .
./mask.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/mask.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/mask.F90 .
./ice_dyn.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_dyn.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_dyn.F90 .
./rot.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/rot.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/rot.F90 .
SRC = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_param/ocean_rough.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/mask.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_type.F90                                          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/rot.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_thm.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_bergs.F90                                                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_grid.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_param/ocean_albedo.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_param/ice_albedo.F90                                       \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_model.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_spec.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_dyn.F90
OBJ = ocean_rough.o mask.o ice_type.o rot.o ice_thm.o ice_bergs.o ice_grid.o ocean_albedo.o ice_albedo.o ice_model.o ice_spec.o ice_dyn.o
OFF = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_param/ocean_rough.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_param/ice_albedo.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_bergs.F90                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_param/ocean_albedo.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_thm.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_spec.F90                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_model.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_type.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_grid.F90                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/mask.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/ice_dyn.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ice_sis/rot.F90
clean: neat
	-rm -f .lib_ice.a.cppdefs $(OBJ) lib_ice.a
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
lib_ice.a: $(OBJ)
	$(AR) $(ARFLAGS) lib_ice.a $(OBJ)
