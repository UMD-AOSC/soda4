# Makefile created by mkmf $Id: mkmf,v 1.1.2.1 2013/12/18 17:47:54 Niki.Zadeh Exp $ 

SRCROOT = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/

CPPDEFS = -Duse_netCDF -Duse_libMPI -DUSE_OCEAN_BGC -DENABLE_ODA -DSPMD -DLAND_BND_TRACERS -DENABLE_GDS -Duse_netCDF4

OTHERFLAGS = -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_FMS -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_ocean -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/exec_test/derecho/lib_version/

include /glade/cheyenne/scratch/chepurin/SODA4_TEST/bin/mkmf.template.derecho


.DEFAULT:
	-echo $@ does not exist.
all: lib_ocean.a
FMS_coupler_util.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/FMS_coupler_util.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/FMS_coupler_util.F90
FMS_ocmip2_co2calc.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/FMS_ocmip2_co2calc.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/FMS_ocmip2_co2calc.F90
adv_center.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/adv_center.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/adv_center.F90
algebraiclength.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/algebraiclength.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/algebraiclength.F90
alpha_mnb.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/alpha_mnb.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                    \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/alpha_mnb.F90
cmue_a.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_a.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_a.F90
cmue_b.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_b.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_b.F90
cmue_c.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_c.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_c.F90
cmue_d.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_d.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_d.F90
cmue_ma.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_ma.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_ma.F90
cmue_rf.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_rf.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_rf.F90
cmue_sg.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_sg.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_sg.F90
compute_cpsi3.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/compute_cpsi3.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/compute_cpsi3.F90
compute_rist.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/compute_rist.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/compute_rist.F90
convert_fluxes.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/convert_fluxes.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h eqstate.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/convert_fluxes.F90
diff_center.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/diff_center.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h util.o tridiagonal_gotm.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/diff_center.F90
diff_face.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/diff_face.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h util.o tridiagonal_gotm.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/diff_face.F90
dissipationeq.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/dissipationeq.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/dissipationeq.F90
epsbalgebraic.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/epsbalgebraic.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/epsbalgebraic.F90
eqstate.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/eqstate.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/eqstate.F90
fk_craig.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/fk_craig.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/fk_craig.F90
generic_BLING.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_BLING.F90 generic_tracer_utils.o FMS_ocmip2_co2calc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_BLING.F90
generic_CFC.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_CFC.F90 generic_tracer_utils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_CFC.F90
generic_COBALT.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_COBALT.F90 generic_tracer_utils.o FMS_ocmip2_co2calc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_COBALT.F90
generic_ERGOM.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_ERGOM.F90 generic_tracer_utils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_ERGOM.F90
generic_TOPAZ.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_TOPAZ.F90 generic_tracer_utils.o FMS_ocmip2_co2calc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_TOPAZ.F90
generic_miniBLING.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_miniBLING.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h generic_tracer_utils.o FMS_ocmip2_co2calc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_miniBLING.F90
generic_tracer.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_tracer.F90 generic_tracer_utils.o generic_CFC.o generic_ERGOM.o generic_TOPAZ.o generic_BLING.o generic_miniBLING.o generic_COBALT.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_tracer.F90
generic_tracer_utils.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_tracer_utils.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h FMS_coupler_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_tracer_utils.F90
genericeq.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/genericeq.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                    \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/genericeq.F90
godas.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/godas.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_types.o ocean_parameters.o godas_types.o godas_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/godas.F90
godas_data.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/godas_data.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/godas_data.F90
godas_types.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/godas_types.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/godas_types.F90
gotm_lib_version.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/gotm_lib_version.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/gotm_lib_version.F90
gridinterpol.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/gridinterpol.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                    \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/gridinterpol.F90
internal_wave.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/internal_wave.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/internal_wave.F90
ispralength.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/ispralength.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/ispralength.F90
kbalgebraic.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/kbalgebraic.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/kbalgebraic.F90
kbeq.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/kbeq.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/kbeq.F90
kpp.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/kpp.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/kpp.F90
lagrange.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/lagrange.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/lagrange.F90
lengthscaleeq.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/lengthscaleeq.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/lengthscaleeq.F90
ocean_adv_vel_diag.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_adv_vel_diag.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_operators.o ocean_parameters.o ocean_types.o \
ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_adv_vel_diag.F90
ocean_advection_velocity.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_advection_velocity.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_operators.o ocean_parameters.o   \
ocean_types.o ocean_workspace.o ocean_obc.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_advection_velocity.F90
ocean_age_tracer.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_age_tracer.F90 ocean_tpm_util.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_age_tracer.F90
ocean_barotropic.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_barotropic.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_bih_friction.o ocean_domains.o ocean_lap_friction.o              \
ocean_obc_barotrop.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_barotropic.F90
ocean_bbc.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_bbc.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_wave_type.o \
ocean_wave.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_bbc.F90
ocean_bgc_restore.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_bgc_restore.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h ocean_tpm_util.o ocean_types.o ocmip2_co2calc.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_bgc_restore.F90
ocean_bih_friction.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bih_friction.F90 ocean_bihcgrid_friction.o ocean_bihcst_friction.o ocean_bihgen_friction.o ocean_domains.o ocean_operators.o ocean_parameters.o           \
ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bih_friction.F90
ocean_bih_tracer.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bih_tracer.F90 ocean_domains.o ocean_operators.o ocean_obc.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bih_tracer.F90
ocean_bihcgrid_friction.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bihcgrid_friction.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_obc.o ocean_operators.o   \
ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bihcgrid_friction.F90
ocean_bihcst_friction.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bihcst_friction.F90 ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bihcst_friction.F90
ocean_bihgen_friction.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bihgen_friction.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_obc.o ocean_operators.o       \
ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bihgen_friction.F90
ocean_blob.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob.F90 ocean_blob_diag.o ocean_blob_dynamic_bottom.o ocean_blob_dynamic_free.o ocean_blob_static_free.o ocean_blob_static_bottom.o ocean_blob_util.o ocean_density.o    \
ocean_domains.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob.F90
ocean_blob_diag.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_diag.F90 ocean_blob_util.o ocean_parameters.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_diag.F90
ocean_blob_dynamic_bottom.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_dynamic_bottom.F90 ocean_blob_util.o ocean_density.o ocean_parameters.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_dynamic_bottom.F90
ocean_blob_dynamic_free.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_dynamic_free.F90 ocean_blob_util.o ocean_density.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_dynamic_free.F90
ocean_blob_static_bottom.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_static_bottom.F90 ocean_blob_util.o ocean_density.o ocean_parameters.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_static_bottom.F90
ocean_blob_static_free.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_static_free.F90 ocean_blob_util.o ocean_density.o ocean_types.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_static_free.F90
ocean_blob_util.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_util.F90 ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_util.F90
ocean_convect.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_convect.F90 ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_convect.F90
ocean_coriolis.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_coriolis.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_parameters.o ocean_types.o ocean_util.o              \
ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_coriolis.F90
ocean_density.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_density.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h                                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_operators.o ocean_parameters.o ocean_pressure.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_density.F90
ocean_diagnostics.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_diagnostics.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_adv_vel_diag.o ocean_domains.o ocean_tracer_diag.o             \
ocean_types.o ocean_velocity_diag.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_diagnostics.F90
ocean_domains.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_domains.F90 ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_domains.F90
ocean_drifters.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_drifters.F90 ocean_domains.o ocean_parameters.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_drifters.F90
ocean_form_drag.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_form_drag.F90 ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_form_drag.F90
ocean_frazil.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_frazil.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_parameters.o ocean_tpm_util.o ocean_types.o           \
ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_frazil.F90
ocean_generic_tracer.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_generic_tracer.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h ocean_tpm_util.o ocean_types.o ocean_parameters.o                \
generic_tracer.o generic_tracer_utils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_generic_tracer.F90
ocean_grids.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_grids.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_grids.F90
ocean_ibgc.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_ibgc.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h ocean_tpm_util.o ocean_types.o ocmip2_co2calc.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_ibgc.F90
ocean_increment_eta.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_increment_eta.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_parameters.o ocean_types.o        \
ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_increment_eta.F90
ocean_increment_tracer.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_increment_tracer.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_parameters.o ocean_types.o  \
ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_increment_tracer.F90
ocean_increment_velocity.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_increment_velocity.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_parameters.o            \
ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_increment_velocity.F90
ocean_lap_friction.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lap_friction.F90 ocean_domains.o ocean_lapcst_friction.o ocean_lapgen_friction.o ocean_lapcgrid_friction.o ocean_operators.o ocean_parameters.o           \
ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lap_friction.F90
ocean_lap_tracer.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lap_tracer.F90 ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lap_tracer.F90
ocean_lapcgrid_friction.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lapcgrid_friction.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_obc.o ocean_operators.o   \
ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lapcgrid_friction.F90
ocean_lapcst_friction.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lapcst_friction.F90 ocean_domains.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lapcst_friction.F90
ocean_lapgen_friction.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lapgen_friction.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_obc.o ocean_operators.o       \
ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lapgen_friction.F90
ocean_mixdownslope.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_mixdownslope.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_density.o ocean_domains.o ocean_parameters.o        \
ocean_types.o ocean_util.o ocean_tracer_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_mixdownslope.F90
ocean_model.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_model.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_advection_velocity.o ocean_barotropic.o ocean_bbc.o ocean_bih_friction.o   \
ocean_blob.o ocean_convect.o ocean_coriolis.o ocean_density.o ocean_diagnostics.o ocean_domains.o ocean_form_drag.o ocean_grids.o ocean_increment_eta.o ocean_increment_tracer.o ocean_increment_velocity.o ocean_lap_tracer.o ocean_bih_tracer.o              \
ocean_lap_friction.o ocean_mixdownslope.o ocean_momentum_source.o ocean_nphysics.o ocean_nphysics_new.o ocean_obc.o ocean_operators.o ocean_overexchange.o ocean_overflow.o ocean_overflow_OFP.o ocean_passive.o ocean_pressure.o ocean_rivermix.o             \
ocean_riverspread.o ocean_parameters.o ocean_sbc.o ocean_shortwave.o ocean_sigma_transport.o ocean_sponges_eta.o ocean_sponges_tracer.o ocean_sponges_velocity.o ocean_submesoscale.o ocean_tempsalt.o ocean_thickness.o ocean_topog.o ocean_tracer_advect.o   \
ocean_tracer.o ocean_tracer_util.o ocean_tpm.o ocean_types.o ocean_util.o ocean_velocity_advect.o ocean_velocity_diag.o ocean_velocity.o ocean_vert_mix.o ocean_vert_gotm.o ocean_workspace.o ocean_xlandinsert.o ocean_xlandmix.o ocean_drifters.o            \
ocean_wave_type.o ocean_wave.o oda_driver.o godas_types.o godas.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_model.F90
ocean_momentum_source.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_momentum_source.F90 ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_momentum_source.F90
ocean_nphysics.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_nphysics_util.o ocean_nphysicsA.o           \
ocean_nphysicsB.o ocean_nphysicsC.o ocean_parameters.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics.F90
ocean_nphysicsA.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysicsA.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_nphysics_util.o ocean_operators.o         \
ocean_parameters.o ocean_sigma_transport.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysicsA.F90
ocean_nphysicsB.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysicsB.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_nphysics_util.o ocean_operators.o         \
ocean_parameters.o ocean_sigma_transport.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysicsB.F90
ocean_nphysicsC.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysicsC.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_density.o ocean_domains.o ocean_nphysics_util.o           \
ocean_operators.o ocean_parameters.o ocean_tracer_diag.o ocean_types.o ocean_util.o ocean_tracer_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysicsC.F90
ocean_nphysics_diff.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_diff.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_nphysics_util_new.o               \
ocean_operators.o ocean_parameters.o ocean_tracer_diag.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_diff.F90
ocean_nphysics_flux.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_flux.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_nphysics_util_new.o               \
ocean_operators.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_flux.F90
ocean_nphysics_new.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_new.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_nphysics_util_new.o                 \
ocean_nphysics_flux.o ocean_nphysics_tensor.o ocean_nphysics_skew.o ocean_nphysics_diff.o ocean_operators.o ocean_parameters.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_new.F90
ocean_nphysics_skew.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_skew.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_nphysics_util_new.o               \
ocean_parameters.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_skew.F90
ocean_nphysics_tensor.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_tensor.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_nphysics_skew.o               \
ocean_sigma_transport.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_tensor.F90
ocean_nphysics_util.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_util.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_density.o ocean_operators.o       \
ocean_parameters.o ocean_tracer_diag.o ocean_types.o ocean_util.o ocean_tracer_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_util.F90
ocean_nphysics_util_new.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_util_new.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_parameters.o              \
ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_util_new.F90
ocean_obc.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_obc.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h ocean_util.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_obc_barotrop.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_obc.F90
ocean_obc_barotrop.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_obc_barotrop.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h ocean_util.o ocean_domains.o ocean_parameters.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_obc_barotrop.F90
ocean_operators.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_operators.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_parameters.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_operators.F90
ocean_overexchange.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_overexchange.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_density.o ocean_domains.o ocean_parameters.o        \
ocean_types.o ocean_util.o ocean_tracer_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_overexchange.F90
ocean_overflow.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_overflow.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o  \
ocean_util.o ocean_tracer_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_overflow.F90
ocean_overflow_OFP.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_overflow_OFP.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_density.o ocean_domains.o ocean_parameters.o        \
ocean_types.o ocean_util.o ocean_tracer_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_overflow_OFP.F90
ocean_parameters.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_parameters.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_parameters.F90
ocean_passive.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_passive.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_tpm_util.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_passive.F90
ocean_pert_co2.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_pert_co2.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h ocean_tpm_util.o ocean_types.o ocmip2_co2calc.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_pert_co2.F90
ocean_po4_pre.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_po4_pre.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h ocean_tpm_util.o ocean_types.o ocean_tracer_diag.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_po4_pre.F90
ocean_pressure.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_pressure.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_operators.o ocean_parameters.o ocean_types.o         \
ocean_util.o ocean_workspace.o ocean_obc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_pressure.F90
ocean_residency.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency.F90 ocean_tpm_util.o ocean_types.o ocean_residency_meta.o ocean_residency_ml.o ocean_residency_range.o ocean_residency_integrand.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency.F90
ocean_residency_integrand.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_integrand.F90 ocean_residency_meta.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_integrand.F90
ocean_residency_meta.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_meta.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_meta.F90
ocean_residency_ml.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_ml.F90 ocean_residency_meta.o ocean_types.o ocean_tracer_diag.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_ml.F90
ocean_residency_range.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_range.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h ocean_residency_meta.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_range.F90
ocean_rivermix.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_rivermix.F90 ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_rivermix.F90
ocean_riverspread.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_riverspread.F90 ocean_domains.o ocean_operators.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_riverspread.F90
ocean_sbc.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_sbc.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h        \
ocean_domains.o ocean_parameters.o ocean_riverspread.o ocean_tempsalt.o ocean_tpm.o ocean_types.o ocean_workspace.o ocean_util.o ocean_tracer_util.o ocean_operators.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_sbc.F90
ocean_shortwave.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h                                                                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_parameters.o ocean_shortwave_csiro.o ocean_shortwave_jerlov.o ocean_shortwave_gfdl.o ocean_types.o ocean_tpm_util.o ocean_workspace.o ocean_util.o        \
ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave.F90
ocean_shortwave_csiro.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave_csiro.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave_csiro.F90
ocean_shortwave_gfdl.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave_gfdl.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_parameters.o ocean_types.o      \
ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave_gfdl.F90
ocean_shortwave_jerlov.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave_jerlov.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_parameters.o ocean_types.o  \
ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave_jerlov.F90
ocean_sigma_transport.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_sigma_transport.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_density.o ocean_operators.o   \
ocean_parameters.o ocean_types.o ocean_util.o ocean_tracer_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_sigma_transport.F90
ocean_sponges_eta.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_sponges_eta.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_parameters.o ocean_types.o            \
ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_sponges_eta.F90
ocean_sponges_tracer.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_sponges_tracer.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_parameters.o ocean_types.o      \
ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_sponges_tracer.F90
ocean_sponges_velocity.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_sponges_velocity.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_parameters.o ocean_types.o  \
ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_sponges_velocity.F90
ocean_submesoscale.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_submesoscale.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_operators.o ocean_parameters.o      \
ocean_tracer_diag.o ocean_types.o ocean_util.o ocean_tracer_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_submesoscale.F90
ocean_tempsalt.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tempsalt.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_parameters.o ocean_tpm_util.o ocean_types.o       \
ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tempsalt.F90
ocean_thickness.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_thickness.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_grids.o ocean_parameters.o ocean_tracer_util.o     \
ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_thickness.F90
ocean_topog.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_topog.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_parameters.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_topog.F90
ocean_tpm.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tpm.F90 ocean_types.o ocean_age_tracer.o ocean_residency.o ocean_pert_co2.o ocmip2_abiotic.o ocmip2_cfc.o ocmip2_biotic.o ocean_bgc_restore.o ocmip2_he.o                \
ocean_po4_pre.o ocean_ibgc.o ocean_generic_tracer.o ocean_frazil.o ocean_tempsalt.o ocean_passive.o transport_matrix.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tpm.F90
ocean_tpm_util.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tpm_util.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tpm_util.F90
ocean_tracer.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tracer.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h transport_matrix.o ocean_blob.o ocean_convect.o ocean_domains.o             \
ocean_density.o ocean_frazil.o ocean_bih_tracer.o ocean_generic_tracer.o ocean_lap_tracer.o ocean_obc.o ocean_parameters.o ocean_passive.o ocean_shortwave.o ocean_tempsalt.o ocean_tpm.o ocean_tpm_util.o ocean_tracer_advect.o ocean_tracer_diag.o           \
ocean_tracer_util.o ocean_thickness.o ocean_types.o ocean_util.o ocean_vert_mix.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tracer.F90
ocean_tracer_advect.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tracer_advect.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_obc.o ocean_parameters.o ocean_topog.o  \
ocean_tracer_util.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tracer_advect.F90
ocean_tracer_diag.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_tracer_diag.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_density.o ocean_domains.o ocean_obc.o ocean_operators.o        \
ocean_parameters.o ocean_tracer_util.o ocean_types.o ocean_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_tracer_diag.F90
ocean_tracer_util.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_tracer_util.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_parameters.o ocean_types.o ocean_util.o        \
ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_tracer_util.F90
ocean_types.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_types.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_types.F90
ocean_util.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_util.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_parameters.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_util.F90
ocean_velocity.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_velocity.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_bih_friction.o ocean_coriolis.o ocean_domains.o ocean_form_drag.o    \
ocean_lap_friction.o ocean_momentum_source.o ocean_obc.o ocean_operators.o ocean_parameters.o ocean_pressure.o ocean_types.o ocean_util.o ocean_velocity_advect.o ocean_velocity_diag.o ocean_vert_mix.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_velocity.F90
ocean_velocity_advect.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_velocity_advect.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_obc.o ocean_operators.o                \
ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_velocity_advect.F90
ocean_velocity_diag.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_velocity_diag.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_bih_friction.o ocean_coriolis.o ocean_domains.o ocean_form_drag.o ocean_lap_friction.o ocean_momentum_source.o ocean_operators.o ocean_parameters.o ocean_pressure.o      \
ocean_types.o ocean_util.o ocean_velocity_advect.o ocean_vert_mix.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_velocity_diag.F90
ocean_vert_chen.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_chen.F90 ocean_density.o ocean_domains.o ocean_parameters.o ocean_shortwave_csiro.o ocean_types.o ocean_util.o ocean_vert_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_chen.F90
ocean_vert_const.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_const.F90 ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_const.F90
ocean_vert_gotm.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_gotm.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h tridiagonal_gotm.o turbulence.o ocean_domains.o                \
ocean_parameters.o ocean_types.o ocean_util.o ocean_workspace.o ocean_obc.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_gotm.F90
ocean_vert_kpp_mom4p0.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_kpp_mom4p0.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_density.o ocean_domains.o ocean_parameters.o \
ocean_types.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_kpp_mom4p0.F90
ocean_vert_kpp_mom4p1.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_kpp_mom4p1.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_density.o ocean_domains.o ocean_parameters.o \
ocean_types.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_kpp_mom4p1.F90
ocean_vert_kpp_test.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_kpp_test.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_density.o ocean_domains.o ocean_parameters.o     \
ocean_types.o ocean_vert_util.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_kpp_test.F90
ocean_vert_mix.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_mix.F90 ocean_domains.o ocean_parameters.o ocean_tracer_util.o ocean_types.o ocean_util.o ocean_vert_const.o ocean_vert_chen.o ocean_vert_gotm.o        \
ocean_vert_kpp_test.o ocean_vert_kpp_mom4p0.o ocean_vert_kpp_mom4p1.o ocean_vert_pp.o ocean_vert_tidal.o ocean_vert_util.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_mix.F90
ocean_vert_pp.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_pp.F90 ocean_density.o ocean_domains.o ocean_parameters.o ocean_types.o ocean_vert_util.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_pp.F90
ocean_vert_tidal.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_tidal.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_operators.o ocean_parameters.o         \
ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_tidal.F90
ocean_vert_tidal_test.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_tidal_test.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_operators.o                  \
ocean_parameters.o ocean_types.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_tidal_test.F90
ocean_vert_util.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_util.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_density.o ocean_domains.o ocean_parameters.o             \
ocean_types.o ocean_workspace.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_util.F90
ocean_wave.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_wave/ocean_wave.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_domains.o ocean_types.o ocean_parameters.o ocean_workspace.o ocean_util.o    \
ocean_wave_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_wave/ocean_wave.F90
ocean_wave_type.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_wave/ocean_wave_type.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h                                                                           \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_wave/ocean_wave_type.F90
ocean_workspace.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_workspace.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_workspace.F90
ocean_xlandinsert.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_xlandinsert.F90 ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_xlandinsert.F90
ocean_xlandmix.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_xlandmix.F90 ocean_domains.o ocean_parameters.o ocean_types.o ocean_workspace.o ocean_util.o ocean_tracer_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_xlandmix.F90
ocmip2_abiotic.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_abiotic.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h ocean_tpm_util.o ocean_types.o ocmip2_co2calc.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_abiotic.F90
ocmip2_biotic.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_biotic.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h ocean_tpm_util.o ocean_types.o ocmip2_co2calc.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_biotic.F90
ocmip2_cfc.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_cfc.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h ocean_tpm_util.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_cfc.F90
ocmip2_co2calc.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_co2calc.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_co2calc.F90
ocmip2_he.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_he.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h ocean_tpm_util.o ocean_types.o ocean_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_he.F90
oda_driver.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/oda_driver.F90 ocean_convect.o ocean_domains.o ocean_grids.o ocean_parameters.o ocean_topog.o ocean_types.o ocean_workspace.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/oda_driver.F90
ode_solvers.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/ode_solvers.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/ode_solvers.F90
potentialml.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/potentialml.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/potentialml.F90
production.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/production.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/production.F90
q2over2eq.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/q2over2eq.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                    \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/q2over2eq.F90
r_ratio.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/r_ratio.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/r_ratio.F90
time.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/time.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                                    \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/time.F90
tkealgebraic.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/tkealgebraic.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/tkealgebraic.F90
tkeeq.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/tkeeq.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/tkeeq.F90
transport_matrix.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/transport_matrix.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h ocean_tpm_util.o ocean_types.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/transport_matrix.F90
tridiagonal_gotm.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/tridiagonal_gotm.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/tridiagonal_gotm.F90
turbulence.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/turbulence.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/turbulence.F90
util.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/util.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                                    \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/util.F90
variances.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/variances.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                    \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h turbulence.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/variances.F90
./ocean_blob_dynamic_free.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_dynamic_free.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_dynamic_free.F90 .
./ocean_sponges_velocity.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_sponges_velocity.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_sponges_velocity.F90 .
./generic_CFC.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_CFC.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_CFC.F90 .
./ocean_nphysicsC.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysicsC.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysicsC.F90 .
./lagrange.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/lagrange.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/lagrange.F90 .
./ocean_sbc.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_sbc.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_sbc.F90 .
./ocean_tempsalt.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tempsalt.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tempsalt.F90 .
./ocean_lapcst_friction.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lapcst_friction.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lapcst_friction.F90 .
./lengthscaleeq.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/lengthscaleeq.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/lengthscaleeq.F90 .
./ocean_ibgc.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_ibgc.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_ibgc.F90 .
./kpp.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/kpp.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/kpp.F90 .
./ocean_lapgen_friction.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lapgen_friction.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lapgen_friction.F90 .
./ocean_generic_tracer.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_generic_tracer.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_generic_tracer.F90 .
./ocean_blob_diag.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_diag.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_diag.F90 .
./diff_face.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/diff_face.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/diff_face.F90 .
./ocean_vert_const.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_const.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_const.F90 .
./ocean_nphysics_tensor.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_tensor.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_tensor.F90 .
./ocean_shortwave_gfdl.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave_gfdl.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave_gfdl.F90 .
./ocean_overexchange.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_overexchange.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_overexchange.F90 .
./gridinterpol.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/gridinterpol.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/gridinterpol.F90 .
./ocmip2_abiotic.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_abiotic.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_abiotic.F90 .
./generic_ERGOM.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_ERGOM.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_ERGOM.F90 .
./ocean_sponges_eta.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_sponges_eta.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_sponges_eta.F90 .
./cppdefs.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h .
./ocean_wave.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_wave/ocean_wave.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_wave/ocean_wave.F90 .
./eqstate.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/eqstate.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/eqstate.F90 .
./ocean_passive.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_passive.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_passive.F90 .
./ocean_velocity_advect.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_velocity_advect.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_velocity_advect.F90 .
./generic_miniBLING.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_miniBLING.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_miniBLING.F90 .
./ocean_bihcst_friction.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bihcst_friction.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bihcst_friction.F90 .
./r_ratio.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/r_ratio.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/r_ratio.F90 .
./ocean_nphysics_diff.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_diff.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_diff.F90 .
./ocean_vert_pp.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_pp.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_pp.F90 .
./ocean_pressure.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_pressure.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_pressure.F90 .
./ispralength.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/ispralength.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/ispralength.F90 .
./ocean_bih_friction.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bih_friction.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bih_friction.F90 .
./genericeq.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/genericeq.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/genericeq.F90 .
./ocean_tpm_util.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tpm_util.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tpm_util.F90 .
./ocean_memory.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h .
./ocean_shortwave.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave.F90 .
./generic_COBALT.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_COBALT.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_COBALT.F90 .
./internal_wave.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/internal_wave.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/internal_wave.F90 .
./ocean_increment_eta.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_increment_eta.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_increment_eta.F90 .
./ocean_pert_co2.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_pert_co2.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_pert_co2.F90 .
./ocean_grids.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_grids.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_grids.F90 .
./alpha_mnb.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/alpha_mnb.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/alpha_mnb.F90 .
./version.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h .
./compute_cpsi3.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/compute_cpsi3.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/compute_cpsi3.F90 .
./ocean_model.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_model.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_model.F90 .
./ocean_momentum_source.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_momentum_source.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_momentum_source.F90 .
./ocean_bihgen_friction.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bihgen_friction.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bihgen_friction.F90 .
./ocean_obc.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_obc.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_obc.F90 .
./algebraiclength.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/algebraiclength.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/algebraiclength.F90 .
./ocmip2_he.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_he.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_he.F90 .
./ocean_nphysics_new.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_new.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_new.F90 .
./ocean_shortwave_csiro.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave_csiro.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave_csiro.F90 .
./kbeq.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/kbeq.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/kbeq.F90 .
./ocean_vert_kpp_test.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_kpp_test.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_kpp_test.F90 .
./cmue_a.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_a.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_a.F90 .
./ocean_velocity_diag.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_velocity_diag.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_velocity_diag.F90 .
./ocean_xlandmix.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_xlandmix.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_xlandmix.F90 .
./time.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/time.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/time.F90 .
./ocean_tracer_advect.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tracer_advect.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tracer_advect.F90 .
./ocean_nphysics_flux.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_flux.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_flux.F90 .
./ocean_blob_dynamic_bottom.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_dynamic_bottom.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_dynamic_bottom.F90 .
./ocean_blob_static_bottom.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_static_bottom.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_static_bottom.F90 .
./ocean_density.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_density.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_density.F90 .
./ocean_riverspread.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_riverspread.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_riverspread.F90 .
./fms_platform.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h .
./ocean_blob.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob.F90 .
./ocmip2_co2calc.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_co2calc.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_co2calc.F90 .
./ocean_tracer.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tracer.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tracer.F90 .
./ocean_wave_type.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_wave/ocean_wave_type.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_wave/ocean_wave_type.F90 .
./ocean_topog.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_topog.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_topog.F90 .
./ocean_bihcgrid_friction.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bihcgrid_friction.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bihcgrid_friction.F90 .
./ocean_vert_kpp_mom4p0.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_kpp_mom4p0.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_kpp_mom4p0.F90 .
./ocean_bih_tracer.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bih_tracer.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bih_tracer.F90 .
./ocean_residency_range.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_range.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_range.F90 .
./FMS_ocmip2_co2calc.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/FMS_ocmip2_co2calc.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/FMS_ocmip2_co2calc.F90 .
./ocean_thickness.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_thickness.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_thickness.F90 .
./generic_BLING.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_BLING.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_BLING.F90 .
./ocean_nphysicsA.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysicsA.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysicsA.F90 .
./ocmip2_biotic.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_biotic.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_biotic.F90 .
./ocean_vert_kpp_mom4p1.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_kpp_mom4p1.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_kpp_mom4p1.F90 .
./ocean_submesoscale.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_submesoscale.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_submesoscale.F90 .
./ocean_vert_mix.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_mix.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_mix.F90 .
./ocean_nphysics.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics.F90 .
./godas.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/godas.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/godas.F90 .
./godas_types.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/godas_types.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/godas_types.F90 .
./tkealgebraic.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/tkealgebraic.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/tkealgebraic.F90 .
./cmue_sg.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_sg.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_sg.F90 .
./ocean_lapcgrid_friction.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lapcgrid_friction.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lapcgrid_friction.F90 .
./ocean_rivermix.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_rivermix.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_rivermix.F90 .
./generic_TOPAZ.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_TOPAZ.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_TOPAZ.F90 .
./ocmip2_cfc.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_cfc.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_cfc.F90 .
./ocean_sponges_tracer.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_sponges_tracer.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_sponges_tracer.F90 .
./turbulence.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/turbulence.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/turbulence.F90 .
./diff_center.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/diff_center.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/diff_center.F90 .
./generic_tracer.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_tracer.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_tracer.F90 .
./ocean_overflow.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_overflow.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_overflow.F90 .
./ocean_xlandinsert.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_xlandinsert.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_xlandinsert.F90 .
./ocean_residency_integrand.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_integrand.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_integrand.F90 .
./ocean_sigma_transport.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_sigma_transport.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_sigma_transport.F90 .
./ocean_blob_static_free.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_static_free.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_static_free.F90 .
./ocean_drifters.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_drifters.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_drifters.F90 .
./convert_fluxes.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/convert_fluxes.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/convert_fluxes.F90 .
./ocean_operators.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_operators.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_operators.F90 .
./ocean_lap_tracer.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lap_tracer.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lap_tracer.F90 .
./ocean_adv_vel_diag.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_adv_vel_diag.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_adv_vel_diag.F90 .
./ocean_increment_velocity.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_increment_velocity.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_increment_velocity.F90 .
./ocean_vert_tidal_test.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_tidal_test.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_tidal_test.F90 .
./cmue_c.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_c.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_c.F90 .
./q2over2eq.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/q2over2eq.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/q2over2eq.F90 .
./ocean_nphysics_skew.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_skew.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_skew.F90 .
./cmue_b.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_b.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_b.F90 .
./ocean_util.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_util.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_util.F90 .
./transport_matrix.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/transport_matrix.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/transport_matrix.F90 .
./ocean_tracer_diag.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_tracer_diag.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_tracer_diag.F90 .
./ocean_tpm.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tpm.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tpm.F90 .
./ocean_nphysicsB.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysicsB.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysicsB.F90 .
./ocean_overflow_OFP.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_overflow_OFP.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_overflow_OFP.F90 .
./ocean_advection_velocity.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_advection_velocity.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_advection_velocity.F90 .
./ocean_nphysics_util_new.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_util_new.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_util_new.F90 .
./ocean_age_tracer.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_age_tracer.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_age_tracer.F90 .
./oda_driver.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/oda_driver.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/oda_driver.F90 .
./adv_center.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/adv_center.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/adv_center.F90 .
./cmue_ma.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_ma.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_ma.F90 .
./ocean_workspace.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_workspace.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_workspace.F90 .
./ocean_vert_chen.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_chen.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_chen.F90 .
./cmue_d.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_d.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_d.F90 .
./ode_solvers.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/ode_solvers.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/ode_solvers.F90 .
./ocean_velocity.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_velocity.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_velocity.F90 .
./cmue_rf.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_rf.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_rf.F90 .
./ocean_convect.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_convect.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_convect.F90 .
./ocean_nphysics_util.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_util.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_util.F90 .
./ocean_bgc_restore.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_bgc_restore.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_bgc_restore.F90 .
./ocean_obc_barotrop.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_obc_barotrop.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_obc_barotrop.F90 .
./ocean_types.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_types.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_types.F90 .
./tridiagonal_gotm.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/tridiagonal_gotm.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/tridiagonal_gotm.F90 .
./ocean_shortwave_jerlov.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave_jerlov.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave_jerlov.F90 .
./ocean_vert_tidal.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_tidal.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_tidal.F90 .
./ocean_coriolis.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_coriolis.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_coriolis.F90 .
./ocean_blob_util.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_util.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_util.F90 .
./ocean_domains.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_domains.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_domains.F90 .
./FMS_coupler_util.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/FMS_coupler_util.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/FMS_coupler_util.F90 .
./ocean_barotropic.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_barotropic.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_barotropic.F90 .
./util.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/util.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/util.F90 .
./epsbalgebraic.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/epsbalgebraic.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/epsbalgebraic.F90 .
./ocean_po4_pre.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_po4_pre.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_po4_pre.F90 .
./fk_craig.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/fk_craig.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/fk_craig.F90 .
./ocean_bbc.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_bbc.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_bbc.F90 .
./generic_tracer_utils.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_tracer_utils.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_tracer_utils.F90 .
./kbalgebraic.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/kbalgebraic.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/kbalgebraic.F90 .
./godas_data.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/godas_data.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/godas_data.F90 .
./ocean_mixdownslope.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_mixdownslope.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_mixdownslope.F90 .
./dissipationeq.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/dissipationeq.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/dissipationeq.F90 .
./compute_rist.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/compute_rist.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/compute_rist.F90 .
./ocean_residency.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency.F90 .
./ocean_frazil.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_frazil.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_frazil.F90 .
./gotm_lib_version.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/gotm_lib_version.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/gotm_lib_version.F90 .
./ocean_parameters.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_parameters.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_parameters.F90 .
./ocean_residency_meta.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_meta.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_meta.F90 .
./ocean_vert_gotm.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_gotm.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_gotm.F90 .
./ocean_lap_friction.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lap_friction.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lap_friction.F90 .
./ocean_residency_ml.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_ml.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_ml.F90 .
./ocean_form_drag.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_form_drag.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_form_drag.F90 .
./production.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/production.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/production.F90 .
./ocean_diagnostics.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_diagnostics.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_diagnostics.F90 .
./variances.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/variances.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/variances.F90 .
./ocean_vert_util.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_util.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_util.F90 .
./ocean_tracer_util.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_tracer_util.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_tracer_util.F90 .
./ocean_increment_tracer.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_increment_tracer.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_increment_tracer.F90 .
./potentialml.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/potentialml.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/potentialml.F90 .
./tkeeq.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/tkeeq.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/tkeeq.F90 .
SRC = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_grids.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_ma.F90                                                                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_ibgc.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_diagnostics.F90                                                                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/kbeq.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/ispralength.F90                                                             \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/dissipationeq.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lapcst_friction.F90                                                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_skew.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_form_drag.F90                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_mixdownslope.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_util.F90                                                                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lapcgrid_friction.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/gridinterpol.F90                                                           \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_CFC.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lapgen_friction.F90                                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/ode_solvers.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_BLING.F90                                                                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_increment_tracer.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bihcgrid_friction.F90                                                       \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_workspace.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_util.F90                                                                           \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/epsbalgebraic.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_new.F90                                                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_velocity.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/diff_face.F90                                                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_barotropic.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave_gfdl.F90                                                                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_range.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_obc_barotrop.F90                                                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_adv_vel_diag.F90                                                                             \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_pressure.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_tracer_diag.F90                                                                                       \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_generic_tracer.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_b.F90                                                                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_gotm.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_util_new.F90                                                             \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_wave/ocean_wave.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/fk_craig.F90                                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_frazil.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_pp.F90                                                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/kpp.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_tidal.F90                                                                    \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_pert_co2.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/godas_types.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_sbc.F90                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_xlandinsert.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_tracer.F90                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_increment_eta.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_COBALT.F90                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_rf.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_a.F90                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_xlandmix.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bihgen_friction.F90                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/tkeeq.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_po4_pre.F90                                                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/adv_center.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_c.F90                                                                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_wave/ocean_wave_type.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/gotm_lib_version.F90                                                                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_rivermix.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_overflow.F90                                                                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_convect.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_domains.F90                                                                                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/lagrange.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_kpp_mom4p0.F90                                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/kbalgebraic.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_chen.F90                                                             \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tracer_advect.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tempsalt.F90                                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_meta.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_age_tracer.F90                                                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/tkealgebraic.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_parameters.F90                                                                     \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_bgc_restore.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_flux.F90                                                                          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/FMS_coupler_util.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_operators.F90                                                                          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_sponges_eta.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave_csiro.F90                                                              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_util.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tpm.F90   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_mix.F90                                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/algebraiclength.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysicsA.F90                                                          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/oda_driver.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/compute_rist.F90                                                                           \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/genericeq.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bih_friction.F90                                                             \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bih_tracer.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_tensor.F90                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_d.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_velocity_diag.F90                                                                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_abiotic.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/alpha_mnb.F90                                                                           \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_passive.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_kpp_test.F90                                                                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_dynamic_bottom.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bihcst_friction.F90                                                              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_advection_velocity.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_diag.F90                                                                              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/diff_center.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_dynamic_free.F90                                                                    \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/production.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_increment_velocity.F90                                                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_model.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/FMS_ocmip2_co2calc.F90                                                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/internal_wave.F90                                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_bbc.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_obc.F90                                                                                                    \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_tracer_utils.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_submesoscale.F90                                                          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_biotic.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_sg.F90                                                                              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lap_friction.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_sponges_velocity.F90                                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tracer.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_coriolis.F90                                                                                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_const.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_types.F90                                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave_jerlov.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysicsC.F90                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_miniBLING.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_drifters.F90                                                                          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/potentialml.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_sponges_tracer.F90                                                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_density.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_cfc.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_velocity_advect.F90      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_ERGOM.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/r_ratio.F90                                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_riverspread.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/util.F90                                                                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/tridiagonal_gotm.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/godas.F90                                                                                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_co2calc.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_diff.F90                                                                             \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/turbulence.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_TOPAZ.F90                                                             \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_ml.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lap_tracer.F90                                                                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_static_free.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysicsB.F90                                                                       \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_thickness.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_he.F90                                                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_overflow_OFP.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/godas_data.F90                                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/time.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/transport_matrix.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_util.F90 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/lengthscaleeq.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_static_bottom.F90                                                           \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_integrand.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_momentum_source.F90                                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/convert_fluxes.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/q2over2eq.F90                                                           \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tpm_util.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/eqstate.F90                                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/compute_cpsi3.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_tracer_util.F90                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_kpp_mom4p1.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_overexchange.F90                                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_topog.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_sigma_transport.F90                                                                             \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_tidal_test.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/variances.F90                                                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h       \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h
OBJ = ocean_grids.o cmue_ma.o ocean_ibgc.o ocean_diagnostics.o kbeq.o ispralength.o dissipationeq.o ocean_lapcst_friction.o ocean_nphysics_skew.o ocean_form_drag.o ocean_mixdownslope.o ocean_vert_util.o ocean_lapcgrid_friction.o gridinterpol.o            \
generic_CFC.o ocean_lapgen_friction.o ode_solvers.o generic_BLING.o ocean_increment_tracer.o ocean_bihcgrid_friction.o ocean_workspace.o ocean_nphysics_util.o epsbalgebraic.o ocean_nphysics_new.o ocean_velocity.o diff_face.o ocean_barotropic.o            \
ocean_shortwave_gfdl.o ocean_residency_range.o ocean_obc_barotrop.o ocean_nphysics.o ocean_adv_vel_diag.o ocean_pressure.o ocean_tracer_diag.o ocean_generic_tracer.o cmue_b.o ocean_vert_gotm.o ocean_nphysics_util_new.o ocean_wave.o fk_craig.o             \
ocean_frazil.o ocean_vert_pp.o kpp.o ocean_vert_tidal.o ocean_pert_co2.o godas_types.o ocean_sbc.o ocean_xlandinsert.o generic_tracer.o ocean_increment_eta.o generic_COBALT.o cmue_rf.o cmue_a.o ocean_xlandmix.o ocean_bihgen_friction.o tkeeq.o             \
ocean_po4_pre.o adv_center.o cmue_c.o ocean_wave_type.o gotm_lib_version.o ocean_rivermix.o ocean_overflow.o ocean_convect.o ocean_domains.o lagrange.o ocean_vert_kpp_mom4p0.o kbalgebraic.o ocean_vert_chen.o ocean_tracer_advect.o ocean_tempsalt.o         \
ocean_residency_meta.o ocean_age_tracer.o tkealgebraic.o ocean_parameters.o ocean_bgc_restore.o ocean_nphysics_flux.o FMS_coupler_util.o ocean_operators.o ocean_sponges_eta.o ocean_shortwave_csiro.o ocean_blob_util.o ocean_residency.o ocean_tpm.o         \
ocean_blob.o ocean_vert_mix.o algebraiclength.o ocean_nphysicsA.o oda_driver.o compute_rist.o genericeq.o ocean_bih_friction.o ocean_bih_tracer.o ocean_nphysics_tensor.o cmue_d.o ocean_velocity_diag.o ocmip2_abiotic.o alpha_mnb.o ocean_passive.o          \
ocean_vert_kpp_test.o ocean_blob_dynamic_bottom.o ocean_bihcst_friction.o ocean_advection_velocity.o ocean_blob_diag.o diff_center.o ocean_blob_dynamic_free.o production.o ocean_increment_velocity.o ocean_model.o FMS_ocmip2_co2calc.o ocean_shortwave.o    \
internal_wave.o ocean_bbc.o ocean_obc.o generic_tracer_utils.o ocean_submesoscale.o ocmip2_biotic.o cmue_sg.o ocean_lap_friction.o ocean_sponges_velocity.o ocean_tracer.o ocean_coriolis.o ocean_vert_const.o ocean_types.o ocean_shortwave_jerlov.o          \
ocean_nphysicsC.o generic_miniBLING.o ocean_drifters.o potentialml.o ocean_sponges_tracer.o ocean_density.o ocmip2_cfc.o ocean_velocity_advect.o generic_ERGOM.o r_ratio.o ocean_riverspread.o util.o tridiagonal_gotm.o godas.o ocmip2_co2calc.o              \
ocean_nphysics_diff.o turbulence.o generic_TOPAZ.o ocean_residency_ml.o ocean_lap_tracer.o ocean_blob_static_free.o ocean_nphysicsB.o ocean_thickness.o ocmip2_he.o ocean_overflow_OFP.o godas_data.o time.o transport_matrix.o ocean_util.o lengthscaleeq.o   \
ocean_blob_static_bottom.o ocean_residency_integrand.o ocean_momentum_source.o convert_fluxes.o q2over2eq.o ocean_tpm_util.o eqstate.o compute_cpsi3.o ocean_tracer_util.o ocean_vert_kpp_mom4p1.o ocean_overexchange.o ocean_topog.o ocean_sigma_transport.o  \
ocean_vert_tidal_test.o variances.o
OFF = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_dynamic_free.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_sponges_velocity.F90                                                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_CFC.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysicsC.F90                                                                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/lagrange.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_sbc.F90                                                                                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tempsalt.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lapcst_friction.F90                                                                       \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/lengthscaleeq.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_ibgc.F90                                                                           \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/kpp.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lapgen_friction.F90                                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_generic_tracer.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_diag.F90                                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/diff_face.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_const.F90                                                                    \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_tensor.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave_gfdl.F90                                                           \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_overexchange.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/gridinterpol.F90                                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_abiotic.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_ERGOM.F90                                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_sponges_eta.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/cppdefs.h                                                                     \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_wave/ocean_wave.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/eqstate.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_passive.F90 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_velocity_advect.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_miniBLING.F90                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bihcst_friction.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/r_ratio.F90                                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_diff.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_pp.F90                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_pressure.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/ispralength.F90                                                                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bih_friction.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/genericeq.F90                                                             \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tpm_util.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_memory.h                                                                                           \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_COBALT.F90                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/internal_wave.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_increment_eta.F90                                                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_pert_co2.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_grids.F90                                                                                              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/alpha_mnb.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/include/version.h                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/compute_cpsi3.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_model.F90                                                                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_momentum_source.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bihgen_friction.F90                                                          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_obc.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/algebraiclength.F90                                                                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_he.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_new.F90                                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave_csiro.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/kbeq.F90                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_kpp_test.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_a.F90                                                              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_velocity_diag.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_xlandmix.F90                                                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/time.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tracer_advect.F90                                                                             \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_flux.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_dynamic_bottom.F90                                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_static_bottom.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_density.F90                                                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_riverspread.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob.F90   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_co2calc.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tracer.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_wave/ocean_wave_type.F90      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_topog.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bihcgrid_friction.F90                                                                           \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_kpp_mom4p0.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_bih_tracer.F90                                                              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_range.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/FMS_ocmip2_co2calc.F90                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_thickness.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_BLING.F90                                                                             \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysicsA.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_biotic.F90                                                                                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_kpp_mom4p1.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_submesoscale.F90                                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_mix.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics.F90                                                                       \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/godas.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/godas_types.F90                                                                                                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/tkealgebraic.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_sg.F90                                                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lapcgrid_friction.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_rivermix.F90                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_TOPAZ.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocmip2_cfc.F90                                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_sponges_tracer.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/turbulence.F90                                                          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/diff_center.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_tracer.F90                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_overflow.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_xlandinsert.F90                                                                     \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_integrand.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_sigma_transport.F90                                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_static_free.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_drifters.F90                                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/convert_fluxes.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_operators.F90                                                                          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lap_tracer.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_adv_vel_diag.F90                                                                           \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_increment_velocity.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_tidal_test.F90                                                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_c.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/q2over2eq.F90                                                             \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_skew.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_b.F90                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_util.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/transport_matrix.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_tracer_diag.F90   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_tpm.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysicsB.F90                                                                                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_overflow_OFP.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_advection_velocity.F90                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_util_new.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_age_tracer.F90                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/oda_driver.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/adv_center.F90                                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_ma.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_workspace.F90                                                                           \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_chen.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_d.F90                                                                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/ode_solvers.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_velocity.F90                                                                              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/cmue_rf.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_convect.F90                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/neutral/ocean_nphysics_util.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_bgc_restore.F90                                                                          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_obc_barotrop.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_types.F90                                                                                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/tridiagonal_gotm.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_shortwave_jerlov.F90                                                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_tidal.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_coriolis.F90                                                                              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_blobs/ocean_blob_util.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_domains.F90                                                                                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/FMS_coupler_util.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_barotropic.F90                                                                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/util/util.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/epsbalgebraic.F90                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_bgc/ocean_po4_pre.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/fk_craig.F90                                                                             \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_bbc.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/ocean_shared/generic_tracers/generic_tracer_utils.F90                                                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/kbalgebraic.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/godas_data.F90                                                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_mixdownslope.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/dissipationeq.F90                                                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/compute_rist.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency.F90                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_frazil.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/gotm_lib_version.F90                                                                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_core/ocean_parameters.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_meta.F90                                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_gotm.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/lateral/ocean_lap_friction.F90                                                                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_tracers/ocean_residency_ml.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_form_drag.F90                                                                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/production.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_diagnostics.F90                                                                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/variances.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/vertical/ocean_vert_util.F90                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_diag/ocean_tracer_util.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/sources/ocean_increment_tracer.F90                                                                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/potentialml.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/mom5/ocean_param/gotm-4.0/turbulence/tkeeq.F90
clean: neat
	-rm -f .lib_ocean.a.cppdefs $(OBJ) lib_ocean.a
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
lib_ocean.a: $(OBJ)
	$(AR) $(ARFLAGS) lib_ocean.a $(OBJ)
