# Makefile created by mkmf $Id: mkmf,v 1.1.2.1 2013/12/18 17:47:54 Niki.Zadeh Exp $ 

SRCROOT = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/

CPPDEFS = -Duse_netCDF -Duse_libMPI -DUSE_OCEAN_BGC -DENABLE_ODA -DSPMD -DLAND_BND_TRACERS -DENABLE_GDS -Duse_netCDF4


include /glade/cheyenne/scratch/chepurin/SODA4_TEST/bin/mkmf.template.derecho


.DEFAULT:
	-echo $@ does not exist.
all: lib_FMS.a
MersenneTwister.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/random_numbers/MersenneTwister.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/random_numbers/MersenneTwister.F90
affinity.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/affinity.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/affinity.c
amip_interp.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/amip_interp/amip_interp.F90 time_interp.o time_manager.o get_cal_time.o mpp_io.o horiz_interp.o fms.o fms_io.o constants.o platform.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/amip_interp/amip_interp.F90
astronomy.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/astronomy/astronomy.F90 fms.o time_manager.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/astronomy/astronomy.F90
atmos_ocean_fluxes.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/coupler/atmos_ocean_fluxes.F90 mpp.o coupler_types.o field_manager.o fm_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/coupler/atmos_ocean_fluxes.F90
axis_utils.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/axis_utils/axis_utils.F90 mpp_io.o mpp.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/axis_utils/axis_utils.F90
block_control.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/block_control/block_control.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/block_control/block_control.F90
cloud_interpolator.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/cloud_interpolator.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/cloud_interpolator.F90
column_diagnostics.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/column_diagnostics/column_diagnostics.F90 mpp_io.o fms.o time_manager.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/column_diagnostics/column_diagnostics.F90
constants.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/constants/constants.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/constants/constants.F90
coupler_types.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/coupler/coupler_types.F90 field_manager.o mpp.o fm_util.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/coupler/coupler_types.F90
create_xgrid.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/create_xgrid.c /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic_util.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/create_xgrid.h          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/create_xgrid.c
data_override.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/data_override/data_override.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h platform.o constants.o mpp_io.o mpp.o horiz_interp.o                     \
time_interp_external.o fms_io.o fms.o axis_utils.o mpp_domains.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/data_override/data_override.F90
diag_axis.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_axis.F90 mpp_domains.o fms.o diag_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_axis.F90
diag_data.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_data.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h time_manager.o mpp_domains.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_data.F90
diag_grid.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_grid.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h constants.o fms.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_grid.F90
diag_manager.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_manager.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h time_manager.o mpp_io.o mpp.o fms.o diag_axis.o diag_util.o diag_data.o     \
diag_table.o diag_output.o diag_grid.o constants.o mpp_domains.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_manager.F90
diag_output.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_output.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h mpp_io.o mpp_domains.o mpp.o diag_axis.o diag_data.o time_manager.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_output.F90
diag_table.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_table.F90 mpp_io.o mpp.o fms.o time_manager.o constants.o diag_data.o diag_util.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_table.F90
diag_util.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_util.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h diag_data.o diag_axis.o diag_output.o diag_grid.o fms.o fms_io.o mpp_domains.o    \
time_manager.o mpp_io.o mpp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_util.F90
drifters.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/fms_switches.h          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_push.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_set_field.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_compute_k.h      \
mpp.o mpp_domains.o drifters_core.o drifters_input.o drifters_io.o drifters_comm.o cloud_interpolator.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters.F90
drifters_comm.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_comm.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h                                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/fms_switches.h mpp.o mpp_domains.o drifters_core.o drifters_input.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_comm.F90
drifters_core.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_core.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_core.F90
drifters_input.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_input.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_input.F90
drifters_io.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_io.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_io.F90
ensemble_manager.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/coupler/ensemble_manager.F90 fms.o mpp.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/coupler/ensemble_manager.F90
fft.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fft/fft.F90 platform.o fms.o fft99.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fft/fft.F90
fft99.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fft/fft99.F90 constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fft/fft99.F90
field_manager.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/field_manager/field_manager.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/field_manager/parse.inc mpp.o mpp_io.o fms.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/field_manager	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/field_manager/field_manager.F90
fm_util.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/field_manager/fm_util.F90 field_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/field_manager/fm_util.F90
fms.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/fms.F90 mpp.o mpp_domains.o mpp_io.o fms_io.o memutils.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/fms.F90
fms_io.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/fms_io.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/read_data_2d.inc                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/read_data_3d.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/read_data_4d.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/write_data.inc mpp_io.o mpp_domains.o mpp.o \
platform.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/fms_io.F90
gaussian_topog.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/topography/gaussian_topog.F90 fms.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/topography/gaussian_topog.F90
get_cal_time.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_manager/get_cal_time.F90 fms.o time_manager.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_manager/get_cal_time.F90
gradient.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/gradient.F90 mpp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/gradient.F90
gradient_c2l.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/gradient_c2l.c /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/constant.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic_util.h              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/gradient_c2l.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/gradient_c2l.c
grid.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/grid.F90 mpp.o constants.o fms.o fms_io.o mosaic.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/grid.F90
horiz_interp.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp.F90 fms.o mpp.o constants.o horiz_interp_type.o horiz_interp_conserve.o horiz_interp_bilinear.o horiz_interp_bicubic.o horiz_interp_spherical.o mpp_io.o      \
mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp.F90
horiz_interp_bicubic.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_bicubic.F90 mpp.o fms.o horiz_interp_type.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_bicubic.F90
horiz_interp_bilinear.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_bilinear.F90 mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_bilinear.F90
horiz_interp_conserve.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_conserve.F90 mpp.o fms.o fms_io.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_conserve.F90
horiz_interp_spherical.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_spherical.F90 mpp.o fms.o constants.o horiz_interp_type.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_spherical.F90
horiz_interp_type.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_type.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_type.F90
interp.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/interp.c /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic_util.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/interp.h                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/create_xgrid.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/interp.c
memuse.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/memutils/memuse.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/memutils/memuse.c
memutils.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/memutils/memutils.F90 mpp.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/memutils/memutils.F90
mosaic.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic.F90 mpp.o mpp_io.o fms_io.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic.F90
mosaic_util.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic_util.c /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic_util.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/constant.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic_util.c
mpp.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/system_clock.h                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_util.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_util_sma.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_util_mpi.inc        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_util_nocomm.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_error_a_a.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_error_a_s.h   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_error_s_a.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_error_s_s.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_comm.inc          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_comm_sma.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_transmit_sma.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_transmit.inc  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_reduce_sma.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_sum_sma.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_sum.inc            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_alltoall_sma.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_comm_mpi.inc                                                                                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_transmit_mpi.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_reduce_mpi.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_sum_mpi.h     \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_alltoall_mpi.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_comm_nocomm.inc                                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_transmit_nocomm.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_reduce_nocomm.h                                                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_sum_nocomm.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_alltoall_nocomm.h                                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_chksum_int.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_chksum_scalar.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_chksum.h     \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_gather.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_scatter.h mpp_parameter.o mpp_data.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp.F90
mpp_data.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_data.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_data_sma.inc          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_data_mpi.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_data_nocomm.inc mpp_parameter.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_data.F90
mpp_domains.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_domains.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h                                                                                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_define_nest_domains.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_util.inc                                                                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_comm.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_define.inc                                                                              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_misc.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_update_domains2D.h                                                                              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_update_domains2D_nonblock.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_update_nonblock.h                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_updateV_nonblock.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_update.h                                                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_updateV.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_check.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_checkV.h       \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_update_nest_domains.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_update_nest.h                                                                           \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_redistribute.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_get_boundary.h                                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_get_boundary.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_group_update.h                                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/group_update_pack.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/group_update_unpack.inc                                                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_reduce.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_global_reduce.h                                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_global_sum.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_global_sum_tl.h                                                                                     \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_global_field.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_global_field.h mpp_parameter.o mpp_data.o mpp.o mpp_memutils.o mpp_pset.o mpp_efp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_domains.F90
mpp_efp.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_efp.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_efp.F90
mpp_io.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_io.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_util.inc               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_misc.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_connect.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_read.inc    \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_read_2Ddecomp.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_read_compressed.h                                                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_read_distributed_ascii.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_read_distributed_ascii.h                                                              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_write.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_write_2Ddecomp.h                                                                                    \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_write_compressed.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_write_unlimited_axis.h                                                                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_write.h mpp_parameter.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_io.F90
mpp_memutils.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_memutils.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_memutils.F90
mpp_parameter.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_parameter.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_parameter.F90
mpp_pset.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_pset.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_pset.F90
mpp_utilities.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_utilities.F90 mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_utilities.F90
nsclock.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/nsclock.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/nsclock.c
oda_core.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/oda_core.F90 fms.o mpp.o mpp_domains.o time_manager.o get_cal_time.o axis_utils.o constants.o oda_types.o write_ocean_data.o mpp_io.o field_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/oda_core.F90
oda_core_ecda.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/oda_core_ecda.F90 fms.o mpp.o mpp_io.o mpp_domains.o mpp_memutils.o time_manager.o get_cal_time.o axis_utils.o horiz_interp_type.o horiz_interp_bilinear.o constants.o       \
oda_types.o xbt_drop_rate_adjust.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/oda_core_ecda.F90
oda_types.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/oda_types.F90 time_manager.o mpp.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/oda_types.F90
platform.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/platform/platform.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/platform/platform.F90
quicksort.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/quicksort.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/quicksort.F90
random_numbers.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/random_numbers/random_numbers.F90 MersenneTwister.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/random_numbers/random_numbers.F90
read_mosaic.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/read_mosaic.c /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/read_mosaic.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/constant.h                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic_util.h
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/read_mosaic.c
sat_vapor_pres.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/sat_vapor_pres/sat_vapor_pres.F90 constants.o fms.o mpp_io.o mpp.o sat_vapor_pres_k.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/sat_vapor_pres/sat_vapor_pres.F90
sat_vapor_pres_k.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/sat_vapor_pres/sat_vapor_pres_k.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/sat_vapor_pres/sat_vapor_pres_k.F90
station_data.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/station_data/station_data.F90 axis_utils.o mpp_io.o fms.o mpp.o mpp_domains.o diag_axis.o diag_output.o diag_manager.o diag_util.o time_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/station_data/station_data.F90
stock_constants.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/exchange/stock_constants.F90 mpp.o time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/exchange/stock_constants.F90
test_fms_io.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/test_fms_io.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h mpp.o mpp_domains.o mpp_io.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/test_fms_io.F90
test_horiz_interp.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/test_horiz_interp.F90 mpp.o mpp_io.o mpp_domains.o horiz_interp.o axis_utils.o fms_io.o fms.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/test_horiz_interp.F90
test_mpp.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h mpp.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp.F90
test_mpp_domains.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp_domains.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h mpp.o mpp_domains.o mpp_memutils.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp_domains.F90
test_mpp_io.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp_io.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h mpp.o mpp_domains.o mpp_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp_io.F90
test_mpp_pset.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp_pset.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h mpp.o mpp_pset.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp_pset.F90
test_xgrid.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/exchange/test_xgrid.F90 mpp.o mpp_domains.o mpp_io.o fms.o fms_io.o constants.o xgrid.o mosaic.o gradient.o ensemble_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/exchange/test_xgrid.F90
threadloc.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/threadloc.c
	$(CC) $(CPPDEFS) $(CPPFLAGS) $(CFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/threadloc.c
time_interp.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_interp/time_interp.F90 time_manager.o fms.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_interp/time_interp.F90
time_interp_external.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_interp/time_interp_external.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h mpp.o mpp_io.o time_manager.o get_cal_time.o mpp_domains.o   \
time_interp.o axis_utils.o fms.o platform.o horiz_interp.o constants.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_interp/time_interp_external.F90
time_manager.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_manager/time_manager.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h constants.o fms.o mpp.o fms_io.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_manager/time_manager.F90
topography.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/topography/topography.F90 gaussian_topog.o horiz_interp.o fms.o fms_io.o constants.o mpp.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/topography/topography.F90
tracer_manager.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/tracer_manager/tracer_manager.F90 mpp.o mpp_io.o fms.o field_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/tracer_manager/tracer_manager.F90
tridiagonal.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/tridiagonal/tridiagonal.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/tridiagonal/tridiagonal.F90
write_ocean_data.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/write_ocean_data.F90 mpp_io.o mpp.o oda_types.o time_manager.o mpp_domains.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/write_ocean_data.F90
xbt_drop_rate_adjust.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/xbt_drop_rate_adjust.f90 oda_types.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) -c	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/xbt_drop_rate_adjust.f90
xgrid.o: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/exchange/xgrid.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h fms.o fms_io.o mpp.o mpp_domains.o mpp_io.o constants.o mosaic.o stock_constants.o gradient.o \
time_manager.o diag_manager.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) -c -I/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include	/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/exchange/xgrid.F90
./mpp_domains_util.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_util.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_util.inc .
./mpp_get_boundary.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_get_boundary.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_get_boundary.h .
./mpp_sum.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_sum.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_sum.inc .
./drifters_set_field.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_set_field.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_set_field.h .
./topography.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/topography/topography.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/topography/topography.F90 .
./horiz_interp_conserve.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_conserve.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_conserve.F90 .
./mpp_utilities.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_utilities.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_utilities.F90 .
./tracer_manager.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/tracer_manager/tracer_manager.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/tracer_manager/tracer_manager.F90 .
./time_interp.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_interp/time_interp.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_interp/time_interp.F90 .
./mosaic_util.c: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic_util.c
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic_util.c .
./horiz_interp_bicubic.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_bicubic.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_bicubic.F90 .
./mpp_read_2Ddecomp.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_read_2Ddecomp.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_read_2Ddecomp.h .
./coupler_types.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/coupler/coupler_types.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/coupler/coupler_types.F90 .
./mpp_util_mpi.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_util_mpi.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_util_mpi.inc .
./mpp_io.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_io.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_io.F90 .
./drifters_push.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_push.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_push.h .
./mpp_do_update_nonblock.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_update_nonblock.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_update_nonblock.h .
./mpp_data_mpi.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_data_mpi.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_data_mpi.inc .
./MersenneTwister.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/random_numbers/MersenneTwister.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/random_numbers/MersenneTwister.F90 .
./mpp_io_read.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_read.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_read.inc .
./diag_grid.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_grid.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_grid.F90 .
./ensemble_manager.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/coupler/ensemble_manager.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/coupler/ensemble_manager.F90 .
./oda_core.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/oda_core.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/oda_core.F90 .
./cloud_interpolator.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/cloud_interpolator.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/cloud_interpolator.F90 .
./mpp_io_misc.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_misc.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_misc.inc .
./mpp_comm_sma.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_comm_sma.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_comm_sma.inc .
./mpp_global_reduce.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_global_reduce.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_global_reduce.h .
./mpp_global_sum_tl.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_global_sum_tl.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_global_sum_tl.h .
./mpp_sum_mpi.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_sum_mpi.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_sum_mpi.h .
./column_diagnostics.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/column_diagnostics/column_diagnostics.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/column_diagnostics/column_diagnostics.F90 .
./mpp_do_global_field.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_global_field.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_global_field.h .
./random_numbers.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/random_numbers/random_numbers.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/random_numbers/random_numbers.F90 .
./write_ocean_data.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/write_ocean_data.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/write_ocean_data.F90 .
./diag_axis.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_axis.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_axis.F90 .
./drifters_comm.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_comm.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_comm.F90 .
./mpp_comm_mpi.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_comm_mpi.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_comm_mpi.inc .
./test_fms_io.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/test_fms_io.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/test_fms_io.F90 .
./xgrid.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/exchange/xgrid.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/exchange/xgrid.F90 .
./mpp_read_compressed.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_read_compressed.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_read_compressed.h .
./mpp_global_field.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_global_field.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_global_field.h .
./group_update_unpack.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/group_update_unpack.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/group_update_unpack.inc .
./constant.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/constant.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/constant.h .
./mpp_read_distributed_ascii.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_read_distributed_ascii.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_read_distributed_ascii.h .
./mpp_chksum.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_chksum.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_chksum.h .
./fms_io.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/fms_io.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/fms_io.F90 .
./mpp_alltoall_nocomm.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_alltoall_nocomm.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_alltoall_nocomm.h .
./mpp_write_compressed.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_write_compressed.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_write_compressed.h .
./stock_constants.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/exchange/stock_constants.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/exchange/stock_constants.F90 .
./atmos_ocean_fluxes.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/coupler/atmos_ocean_fluxes.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/coupler/atmos_ocean_fluxes.F90 .
./mpp_reduce_sma.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_reduce_sma.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_reduce_sma.h .
./mpp_pset.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_pset.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_pset.F90 .
./mpp_domains_reduce.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_reduce.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_reduce.inc .
./mpp_group_update.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_group_update.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_group_update.h .
./horiz_interp_spherical.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_spherical.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_spherical.F90 .
./mpp_domains.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_domains.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_domains.F90 .
./mpp_util_nocomm.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_util_nocomm.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_util_nocomm.inc .
./mpp_data.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_data.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_data.F90 .
./test_mpp_pset.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp_pset.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp_pset.F90 .
./test_mpp_io.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp_io.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp_io.F90 .
./read_mosaic.c: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/read_mosaic.c
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/read_mosaic.c .
./gradient_c2l.c: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/gradient_c2l.c
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/gradient_c2l.c .
./time_interp_external.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_interp/time_interp_external.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_interp/time_interp_external.F90 .
./constants.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/constants/constants.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/constants/constants.F90 .
./mpp_write_2Ddecomp.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_write_2Ddecomp.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_write_2Ddecomp.h .
./memuse.c: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/memutils/memuse.c
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/memutils/memuse.c .
./memutils.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/memutils/memutils.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/memutils/memutils.F90 .
./mpp_write.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_write.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_write.h .
./station_data.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/station_data/station_data.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/station_data/station_data.F90 .
./block_control.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/block_control/block_control.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/block_control/block_control.F90 .
./mpp_io_util.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_util.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_util.inc .
./mpp_do_redistribute.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_redistribute.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_redistribute.h .
./mpp_error_a_s.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_error_a_s.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_error_a_s.h .
./interp.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/interp.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/interp.h .
./mpp_io_write.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_write.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_write.inc .
./diag_data.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_data.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_data.F90 .
./fft99.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fft/fft99.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fft/fft99.F90 .
./mpp_util_sma.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_util_sma.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_util_sma.inc .
./write_data.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/write_data.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/write_data.inc .
./mpp_parameter.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_parameter.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_parameter.F90 .
./mpp_do_get_boundary.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_get_boundary.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_get_boundary.h .
./oda_core_ecda.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/oda_core_ecda.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/oda_core_ecda.F90 .
./mpp_transmit_mpi.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_transmit_mpi.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_transmit_mpi.h .
./system_clock.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/system_clock.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/system_clock.h .
./mpp_efp.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_efp.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_efp.F90 .
./mpp_do_updateV_nonblock.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_updateV_nonblock.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_updateV_nonblock.h .
./diag_manager.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_manager.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_manager.F90 .
./mpp_data_nocomm.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_data_nocomm.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_data_nocomm.inc .
./fms_platform.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h .
./gradient_c2l.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/gradient_c2l.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/gradient_c2l.h .
./create_xgrid.c: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/create_xgrid.c
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/create_xgrid.c .
./grid.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/grid.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/grid.F90 .
./test_xgrid.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/exchange/test_xgrid.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/exchange/test_xgrid.F90 .
./time_manager.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_manager/time_manager.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_manager/time_manager.F90 .
./mpp_scatter.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_scatter.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_scatter.h .
./data_override.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/data_override/data_override.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/data_override/data_override.F90 .
./mpp_do_checkV.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_checkV.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_checkV.h .
./read_data_3d.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/read_data_3d.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/read_data_3d.inc .
./mosaic.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic.F90 .
./test_horiz_interp.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/test_horiz_interp.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/test_horiz_interp.F90 .
./mpp.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp.F90 .
./mpp_io_connect.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_connect.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_connect.inc .
./threadloc.c: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/threadloc.c
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/threadloc.c .
./mpp_error_a_a.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_error_a_a.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_error_a_a.h .
./field_manager.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/field_manager/field_manager.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/field_manager/field_manager.F90 .
./horiz_interp_bilinear.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_bilinear.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_bilinear.F90 .
./fms_switches.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/fms_switches.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/fms_switches.h .
./mpp_comm.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_comm.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_comm.inc .
./mpp_reduce_nocomm.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_reduce_nocomm.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_reduce_nocomm.h .
./mpp_domains_comm.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_comm.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_comm.inc .
./read_data_2d.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/read_data_2d.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/read_data_2d.inc .
./mpp_gather.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_gather.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_gather.h .
./drifters_io.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_io.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_io.F90 .
./get_cal_time.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_manager/get_cal_time.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_manager/get_cal_time.F90 .
./drifters.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters.F90 .
./fms.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/fms.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/fms.F90 .
./mpp_alltoall_sma.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_alltoall_sma.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_alltoall_sma.h .
./fm_util.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/field_manager/fm_util.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/field_manager/fm_util.F90 .
./mpp_alltoall_mpi.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_alltoall_mpi.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_alltoall_mpi.h .
./drifters_core.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_core.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_core.F90 .
./affinity.c: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/affinity.c
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/affinity.c .
./quicksort.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/quicksort.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/quicksort.F90 .
./mpp_error_s_a.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_error_s_a.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_error_s_a.h .
./test_mpp.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp.F90 .
./mpp_sum_sma.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_sum_sma.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_sum_sma.h .
./mpp_transmit_nocomm.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_transmit_nocomm.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_transmit_nocomm.h .
./mpp_sum_nocomm.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_sum_nocomm.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_sum_nocomm.h .
./group_update_pack.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/group_update_pack.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/group_update_pack.inc .
./sat_vapor_pres.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/sat_vapor_pres/sat_vapor_pres.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/sat_vapor_pres/sat_vapor_pres.F90 .
./mpp_chksum_scalar.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_chksum_scalar.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_chksum_scalar.h .
./sat_vapor_pres_k.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/sat_vapor_pres/sat_vapor_pres_k.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/sat_vapor_pres/sat_vapor_pres_k.F90 .
./mpp_reduce_mpi.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_reduce_mpi.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_reduce_mpi.h .
./amip_interp.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/amip_interp/amip_interp.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/amip_interp/amip_interp.F90 .
./mpp_update_domains2D.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_update_domains2D.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_update_domains2D.h .
./platform.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/platform/platform.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/platform/platform.F90 .
./parse.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/field_manager/parse.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/field_manager/parse.inc .
./mpp_define_nest_domains.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_define_nest_domains.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_define_nest_domains.inc .
./gradient.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/gradient.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/gradient.F90 .
./horiz_interp.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp.F90 .
./nsclock.c: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/nsclock.c
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/nsclock.c .
./mpp_global_sum.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_global_sum.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_global_sum.h .
./xbt_drop_rate_adjust.f90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/xbt_drop_rate_adjust.f90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/xbt_drop_rate_adjust.f90 .
./drifters_compute_k.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_compute_k.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_compute_k.h .
./mpp_error_s_s.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_error_s_s.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_error_s_s.h .
./mpp_data_sma.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_data_sma.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_data_sma.inc .
./create_xgrid.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/create_xgrid.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/create_xgrid.h .
./mpp_transmit.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_transmit.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_transmit.inc .
./astronomy.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/astronomy/astronomy.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/astronomy/astronomy.F90 .
./drifters_input.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_input.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_input.F90 .
./mpp_do_update_nest.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_update_nest.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_update_nest.h .
./horiz_interp_type.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_type.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_type.F90 .
./mpp_do_update.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_update.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_update.h .
./mpp_transmit_sma.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_transmit_sma.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_transmit_sma.h .
./read_data_4d.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/read_data_4d.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/read_data_4d.inc .
./mpp_update_nest_domains.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_update_nest_domains.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_update_nest_domains.h .
./diag_util.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_util.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_util.F90 .
./mpp_update_domains2D_nonblock.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_update_domains2D_nonblock.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_update_domains2D_nonblock.h .
./mpp_write_unlimited_axis.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_write_unlimited_axis.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_write_unlimited_axis.h .
./tridiagonal.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/tridiagonal/tridiagonal.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/tridiagonal/tridiagonal.F90 .
./mosaic_util.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic_util.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic_util.h .
./mpp_read_distributed_ascii.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_read_distributed_ascii.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_read_distributed_ascii.inc .
./mpp_comm_nocomm.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_comm_nocomm.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_comm_nocomm.inc .
./mpp_domains_misc.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_misc.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_misc.inc .
./mpp_util.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_util.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_util.inc .
./diag_table.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_table.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_table.F90 .
./test_mpp_domains.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp_domains.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp_domains.F90 .
./interp.c: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/interp.c
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/interp.c .
./axis_utils.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/axis_utils/axis_utils.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/axis_utils/axis_utils.F90 .
./mpp_memutils.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_memutils.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_memutils.F90 .
./mpp_do_check.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_check.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_check.h .
./read_mosaic.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/read_mosaic.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/read_mosaic.h .
./fft.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fft/fft.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fft/fft.F90 .
./gaussian_topog.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/topography/gaussian_topog.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/topography/gaussian_topog.F90 .
./mpp_domains_define.inc: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_define.inc
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_define.inc .
./diag_output.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_output.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_output.F90 .
./mpp_chksum_int.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_chksum_int.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_chksum_int.h .
./oda_types.F90: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/oda_types.F90
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/oda_types.F90 .
./mpp_do_updateV.h: /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_updateV.h
	cp /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_updateV.h .
SRC = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_conserve.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/random_numbers/random_numbers.F90                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/memutils/memutils.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_table.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/station_data/station_data.F90           \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/column_diagnostics/column_diagnostics.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic.F90                                                                                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_type.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/interp.c /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/gradient_c2l.c                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_memutils.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/test_fms_io.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/tracer_manager/tracer_manager.F90                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_interp/time_interp.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/topography/topography.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp_pset.F90               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/oda_core.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/block_control/block_control.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp_domains.F90           \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/data_override/data_override.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/tridiagonal/tridiagonal.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_pset.F90              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_manager/get_cal_time.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/oda_types.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/coupler/coupler_types.F90           \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/xbt_drop_rate_adjust.f90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_comm.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/exchange/xgrid.F90          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/field_manager/field_manager.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_domains.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/topography/gaussian_topog.F90         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/constants/constants.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/axis_utils/axis_utils.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_util.F90              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_data.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_io.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/affinity.c                                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_data.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_output.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/coupler/ensemble_manager.F90      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/exchange/stock_constants.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/platform/platform.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/sat_vapor_pres/sat_vapor_pres.F90      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_manager.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_parameter.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/threadloc.c                       \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/quicksort.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_bilinear.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/gradient.F90         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_io.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/test_horiz_interp.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_bicubic.F90   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fft/fft99.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/field_manager/fm_util.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_core.F90                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_utilities.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/astronomy/astronomy.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fft/fft.F90                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp_io.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/oda_core_ecda.F90                          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_interp/time_interp_external.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/create_xgrid.c                    \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_spherical.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/nsclock.c               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/fms_io.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic_util.c /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_grid.F90                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_axis.F90                                                                                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/sat_vapor_pres/sat_vapor_pres_k.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/fms.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/coupler/atmos_ocean_fluxes.F90            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/read_mosaic.c /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/grid.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_manager/time_manager.F90                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/memutils/memuse.c /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/exchange/test_xgrid.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_efp.F90                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/write_ocean_data.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_input.F90                                                                                       \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/random_numbers/MersenneTwister.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/amip_interp/amip_interp.F90                                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/cloud_interpolator.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_alltoall_sma.h                                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_read.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_gather.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_data_mpi.inc         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_update_nonblock.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_push.h                                                                                    \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_util_mpi.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/read_data_2d.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_comm.inc        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_reduce_nocomm.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_read_2Ddecomp.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_comm.inc  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/fms_switches.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_connect.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_error_a_a.h        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_sum.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_set_field.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_checkV.h         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/read_data_3d.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_scatter.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_get_boundary.h             \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_util.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_alltoall_nocomm.h                                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_define_nest_domains.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_chksum.h                                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_read_distributed_ascii.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/field_manager/parse.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/constant.h        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/group_update_unpack.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_global_field.h                                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_read_compressed.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_update_domains2D.h                                                                             \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_comm_mpi.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_reduce_mpi.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_chksum_scalar.h \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/group_update_pack.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_sum_nocomm.h                                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_transmit_nocomm.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_sum_sma.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_error_s_a.h   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_alltoall_mpi.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_global_field.h                                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_sum_mpi.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_global_sum_tl.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_global_reduce.h \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_misc.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_comm_sma.inc                                                                                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_write_unlimited_axis.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_util.inc                                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_update_domains2D_nonblock.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_update_nest_domains.h                                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_write.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/read_data_4d.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_transmit_sma.h               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_update.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_update_nest.h                                                                                     \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_write_2Ddecomp.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_transmit.inc                                                                                    \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_util_nocomm.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/create_xgrid.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_group_update.h      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_data_sma.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_reduce.inc                                                                                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_compute_k.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_error_s_s.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_reduce_sma.h    \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_global_sum.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_write_compressed.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/gradient_c2l.h     \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_data_nocomm.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_updateV.h       \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_updateV_nonblock.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_chksum_int.h                                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/system_clock.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_transmit_mpi.h                                                                                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_define.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/read_mosaic.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_check.h        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_get_boundary.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_util_sma.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/write_data.inc         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_util.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_misc.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_comm_nocomm.inc \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_write.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/interp.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_read_distributed_ascii.inc   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic_util.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_error_a_s.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_redistribute.h
OBJ = horiz_interp_conserve.o random_numbers.o memutils.o diag_table.o station_data.o column_diagnostics.o mosaic.o horiz_interp_type.o interp.o gradient_c2l.o mpp_memutils.o test_fms_io.o tracer_manager.o time_interp.o topography.o test_mpp_pset.o       \
oda_core.o block_control.o test_mpp_domains.o data_override.o tridiagonal.o mpp_pset.o get_cal_time.o oda_types.o coupler_types.o xbt_drop_rate_adjust.o drifters_comm.o xgrid.o field_manager.o mpp_domains.o gaussian_topog.o constants.o axis_utils.o       \
diag_util.o mpp_data.o drifters_io.o affinity.o diag_data.o diag_output.o ensemble_manager.o stock_constants.o platform.o sat_vapor_pres.o diag_manager.o mpp_parameter.o threadloc.o quicksort.o horiz_interp_bilinear.o gradient.o mpp_io.o                  \
test_horiz_interp.o horiz_interp_bicubic.o fft99.o fm_util.o drifters_core.o mpp_utilities.o astronomy.o fft.o test_mpp.o test_mpp_io.o oda_core_ecda.o mpp.o time_interp_external.o create_xgrid.o horiz_interp_spherical.o drifters.o nsclock.o fms_io.o     \
mosaic_util.o diag_grid.o horiz_interp.o diag_axis.o sat_vapor_pres_k.o fms.o atmos_ocean_fluxes.o read_mosaic.o grid.o time_manager.o memuse.o test_xgrid.o mpp_efp.o write_ocean_data.o drifters_input.o MersenneTwister.o amip_interp.o                     \
cloud_interpolator.o
OFF = /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_util.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_get_boundary.h                                                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_sum.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_set_field.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/topography/topography.F90           \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_conserve.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_utilities.F90                                                                                     \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/tracer_manager/tracer_manager.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_interp/time_interp.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic_util.c        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_bicubic.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_read_2Ddecomp.h                                                                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/coupler/coupler_types.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_util_mpi.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_io.F90                     \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_push.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_update_nonblock.h                                                                                    \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_data_mpi.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/random_numbers/MersenneTwister.F90                                                                                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_read.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_grid.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/coupler/ensemble_manager.F90       \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/oda_core.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/cloud_interpolator.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_misc.inc        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_comm_sma.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_global_reduce.h                                                                                     \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_global_sum_tl.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_sum_mpi.h                                                                                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/column_diagnostics/column_diagnostics.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_global_field.h                                                                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/random_numbers/random_numbers.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/write_ocean_data.F90                                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_axis.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_comm.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_comm_mpi.inc        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/test_fms_io.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/exchange/xgrid.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_read_compressed.h                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_global_field.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/group_update_unpack.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/constant.h      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_read_distributed_ascii.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_chksum.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/fms_io.F90          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_alltoall_nocomm.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_write_compressed.h                                                                             \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/exchange/stock_constants.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/coupler/atmos_ocean_fluxes.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_reduce_sma.h  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_pset.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_reduce.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_group_update.h        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_spherical.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_domains.F90                                                                                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_util_nocomm.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_data.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp_pset.F90                    \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp_io.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/read_mosaic.c /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/gradient_c2l.c                            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_interp/time_interp_external.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/constants/constants.F90                                                                                     \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_write_2Ddecomp.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/memutils/memuse.c /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/memutils/memutils.F90                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_write.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/station_data/station_data.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/block_control/block_control.F90     \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_util.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_redistribute.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_error_a_s.h \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/interp.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_write.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_data.F90                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fft/fft99.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_util_sma.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/write_data.inc                             \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_parameter.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_get_boundary.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/oda_core_ecda.F90       \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_transmit_mpi.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/system_clock.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_efp.F90                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_updateV_nonblock.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_manager.F90                                                                              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_data_nocomm.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/include/fms_platform.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/gradient_c2l.h              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/create_xgrid.c /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/grid.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/exchange/test_xgrid.F90                             \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_manager/time_manager.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_scatter.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/data_override/data_override.F90   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_checkV.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/read_data_3d.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic.F90                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/test_horiz_interp.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_io_connect.inc             \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/threadloc.c /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_error_a_a.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/field_manager/field_manager.F90               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_bilinear.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/fms_switches.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_comm.inc   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_reduce_nocomm.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_comm.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/read_data_2d.inc     \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_gather.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_io.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/time_manager/get_cal_time.F90           \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/fms.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_alltoall_sma.h                          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/field_manager/fm_util.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_alltoall_mpi.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_core.F90       \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/affinity.c /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/quicksort.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_error_s_a.h                         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_sum_sma.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_transmit_nocomm.h              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_sum_nocomm.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/group_update_pack.inc                                                                                   \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/sat_vapor_pres/sat_vapor_pres.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_chksum_scalar.h                                                                                \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/sat_vapor_pres/sat_vapor_pres_k.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_reduce_mpi.h                                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/amip_interp/amip_interp.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_update_domains2D.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/platform/platform.F90      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/field_manager/parse.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_define_nest_domains.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/gradient.F90       \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/nsclock.c /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_global_sum.h                  \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/xbt_drop_rate_adjust.f90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_compute_k.h                                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_error_s_s.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_data_sma.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/create_xgrid.h            \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_transmit.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/astronomy/astronomy.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/drifters/drifters_input.F90          \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_update_nest.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/horiz_interp/horiz_interp_type.F90                                                                              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_update.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_transmit_sma.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fms/read_data_4d.inc           \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_update_nest_domains.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_util.F90                                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_update_domains2D_nonblock.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_write_unlimited_axis.h                                                               \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/tridiagonal/tridiagonal.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/mosaic_util.h                                                                                                 \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_read_distributed_ascii.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_comm_nocomm.inc                                                                       \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_misc.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_util.inc /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_table.F90     \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/test_mpp_domains.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/interp.c /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/axis_utils/axis_utils.F90                        \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/mpp_memutils.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_check.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mosaic/read_mosaic.h                      \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/fft/fft.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/topography/gaussian_topog.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_domains_define.inc              \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/diag_manager/diag_output.F90 /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_chksum_int.h /glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/oda_tools/oda_types.F90         \
/glade/cheyenne/scratch/chepurin/SODA4_TEST/src/shared/mpp/include/mpp_do_updateV.h
clean: neat
	-rm -f .lib_FMS.a.cppdefs $(OBJ) lib_FMS.a
neat:
	-rm -f $(TMPFILES)
localize: $(OFF)
	cp $(OFF) .
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
lib_FMS.a: $(OBJ)
	$(AR) $(ARFLAGS) lib_FMS.a $(OBJ)
