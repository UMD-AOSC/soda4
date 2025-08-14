# Makefile created by mkmf $Id: mkmf,v 1.1.2.1 2013/12/18 17:47:54 Niki.Zadeh Exp $ 

SRCROOT = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/

CPPDEFS = -Duse_netCDF -Duse_libMPI -DUSE_OCEAN_BGC -DENABLE_ODA -DSPMD -DLAND_BND_TRACERS -DENABLE_GDS -Duse_netCDF4

OTHERFLAGS = -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_FMS

include /glade/cheyenne/scratch/chepurin/SODA4_TEST/bin/mkmf.template.derecho


.DEFAULT:
	-echo $@ does not exist.
all: lib_atmos_null.a
atmos_model.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_null/atmos_model.F90 diag_integral.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_null/atmos_model.F90
atmos_tracer_driver.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_null/atmos_tracer_driver.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_null/atmos_tracer_driver.F90
diag_integral.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/diag_integral/diag_integral.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/diag_integral/diag_integral.F90
monin_obukhov.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/monin_obukhov/monin_obukhov.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/monin_obukhov/monin_obukhov.F90
monin_obukhov_kernel.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/monin_obukhov/monin_obukhov_kernel.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h                                                       \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/monin_obukhov/monin_obukhov_interfaces.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/monin_obukhov	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/monin_obukhov/monin_obukhov_kernel.F90
./monin_obukhov_interfaces.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/monin_obukhov/monin_obukhov_interfaces.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/monin_obukhov/monin_obukhov_interfaces.h .
./diag_integral.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/diag_integral/diag_integral.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/diag_integral/diag_integral.F90 .
./atmos_tracer_driver.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_null/atmos_tracer_driver.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_null/atmos_tracer_driver.F90 .
./atmos_model.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_null/atmos_model.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_null/atmos_model.F90 .
./monin_obukhov_kernel.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/monin_obukhov/monin_obukhov_kernel.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/monin_obukhov/monin_obukhov_kernel.F90 .
./monin_obukhov.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/monin_obukhov/monin_obukhov.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/monin_obukhov/monin_obukhov.F90 .
./fms_platform.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h .
SRC = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_null/atmos_model.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/monin_obukhov/monin_obukhov_kernel.F90                                                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/monin_obukhov/monin_obukhov.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_null/atmos_tracer_driver.F90                                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/diag_integral/diag_integral.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h                                                                                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/monin_obukhov/monin_obukhov_interfaces.h
OBJ = atmos_model.o monin_obukhov_kernel.o monin_obukhov.o atmos_tracer_driver.o diag_integral.o
OFF = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/monin_obukhov/monin_obukhov_interfaces.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/diag_integral/diag_integral.F90                                                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_null/atmos_tracer_driver.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_null/atmos_model.F90                                                                                                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/monin_obukhov/monin_obukhov_kernel.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/atmos_param/monin_obukhov/monin_obukhov.F90                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h
clean: neat
	-rm -f .lib_atmos_null.a.cppdefs $(OBJ) lib_atmos_null.a
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
lib_atmos_null.a: $(OBJ)
	$(AR) $(ARFLAGS) lib_atmos_null.a $(OBJ)
