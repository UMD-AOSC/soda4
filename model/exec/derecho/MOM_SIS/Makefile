# Makefile created by mkmf $Id: mkmf,v 1.1.2.1 2013/12/18 17:47:54 Niki.Zadeh Exp $ 

SRCROOT = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/

CPPDEFS = -Duse_netCDF -Duse_libMPI -DUSE_OCEAN_BGC -DENABLE_ODA -DSPMD -DLAND_BND_TRACERS -DENABLE_GDS -Duse_netCDF4

OTHERFLAGS = -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_FMS -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_ocean -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_version/ -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_ice -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_atmos_null -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_land_null

include /glade/cheyenne/scratch/chepurin/SODA4_TEST/bin/mkmf.template.derecho


.DEFAULT:
	-echo $@ does not exist.
all: fms_MOM_SIS.x
coupler_main.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/coupler/coupler_main.F90 flux_exchange.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/coupler/coupler_main.F90
flux_exchange.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/coupler/flux_exchange.F90 surface_flux.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/coupler/flux_exchange.F90
surface_flux.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/coupler/surface_flux.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/coupler/surface_flux.F90
./surface_flux.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/coupler/surface_flux.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/coupler/surface_flux.F90 .
./flux_exchange.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/coupler/flux_exchange.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/coupler/flux_exchange.F90 .
./coupler_main.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/coupler/coupler_main.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/coupler/coupler_main.F90 .
SRC = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/coupler/surface_flux.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/coupler/flux_exchange.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/coupler/coupler_main.F90
OBJ = surface_flux.o flux_exchange.o coupler_main.o
OFF = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/coupler/surface_flux.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/coupler/flux_exchange.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/coupler/coupler_main.F90
clean: neat
	-rm -f .fms_MOM_SIS.x.cppdefs $(OBJ) fms_MOM_SIS.x
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
fms_MOM_SIS.x: $(OBJ) /glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_ocean/lib_ocean.a /glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_ice/lib_ice.a /glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_atmos_null/lib_atmos_null.a /glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_land_null/lib_land_null.a /glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_version/lib_version.a /glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_FMS/lib_FMS.a
	$(LD) $(OBJ) -o fms_MOM_SIS.x /glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_ocean/lib_ocean.a /glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_ice/lib_ice.a /glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_atmos_null/lib_atmos_null.a /glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_land_null/lib_land_null.a /glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_version/lib_version.a /glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_FMS/lib_FMS.a $(LDFLAGS)
