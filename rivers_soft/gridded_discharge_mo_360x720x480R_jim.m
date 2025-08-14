% This routine reads the Dai(2017) monthly river discharge data 1980-2019, filling missing data with the 1980-2017 monthly
% climatology, replaces seven Arctic river discharge data sets with the Great Arctic Rivers data, and adds the Bamber et al (2019) 
% total Greenland discharge. Then the data is mapped onto a 1/2-deg global monthly grid (with zeros where no data is present) and 
% written out in a nc4 file: gridded_discharge_mo_360x720x480R.nc
%
% calls functions: get_rivers, get_Arctic_rivers, and Add_Greenland
%
close all hidden; % clear everything
addpath('/aosc/horse/carton/rivers','-end');
% Get the Dai-Trenberth river discharge data set for 1/1980-12/2017 (extended to 2019, 480 months) including filling missing data with climatological monthly data
%
[fa_full,fa_clim,fa_lon,fa_lat] = get_rivers;
%
[Aa_full,fa_full2] = get_Arctic_rivers(fa_full); % loads Aa_full with Arctic river discharge time series [480x7] and updates fa_full
clear fa_full;
%
[Ba_total_discharge_80_19,Ba_lon,Ba_lat]= Add_Greenland; % loads Ba_total_discharge_80_19 with Greenland discharge time series [752x785x480] (units:m3/s)
% 
%  initialize a 1/2x1/2-deg discharge field spanning [90S-90N, 0E-360E] and [1/1980 to 12/2019]. 
%
gridded_discharg = zeros(480,360,720); 
%
% outer loop to load the 925 river discharge time series
%
for irev = 1:size(fa_full2,2)
% convert discharge locations [lat,lon] to gridpoint space [j,i] by rounding and adding 1 (to guarantee i,j >= 1)
	i = (round(fa_lon(irev,1)*2)+360)+1; % longitudes are given in the range [-180 to 180E] while we shift to a grid where i=1 corresponds to 0E
	j = (round(fa_lat(irev,1)*2)+180)+1;
% inner loop over time
	for itime = 1:480 % itime=1 corresponds to 1/1980, itime=480 corresponds to 12/2019
		gridded_discharg(itime,j,i) = fa_full2(itime,irev);
end;end;
%
% scan file Ba_total_discharge_80_19 through all 752x785 points at each time 1:480 and add discharge to file gridded_discharg
%
for irev = 1:size(Ba_total_discharge_80_19,1);
	for jrev = 1:size(Ba_total_discharge_80_19,2);
% convert discharge locations [lat,lon] to gridpoint space [j,i] by rounding and adding 1 (to guarantee i,j >= 1)
			i = (round(Ba_lon(irev,jrev)*2)+360)+1; % longitudes range: [-132.882 to +43.0501] while we shift to a grid where i=1 corresponds to 0E
			j = (round(Ba_lat(irev,jrev)*2)+180)+1; % latitudes range: [50.4958 to 73.6820] 
% inner loop over time
			for itime = 1:480 % itime=1 corresponds to 1/1980, itime=480 corresponds to 12/2019
				gridded_discharg(itime,j,i) = gridded_discharg(itime,j,i)+Ba_total_discharge_80_19(irev,jrev,itime);
end;end;end;
% Check: contourf(squeeze(mean(gridded_discharg,1)))
% check: mean(sum(gridded_discharg(373:480,:,:), [2 3]),1) 3x104 m3/s for 1980-2019, 3.6x104 m3/2 for 2010-2019
%

for i = 1:720;
	for j=1:360;
		x(i)=(i-361)/2.; %this inverts line #60
		y(j)=(j-181)/2;%this inverts line #60
end;end;
for it = 1:480
	timenum(it)=it;
end;
%
nccreate('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','time','dimensions',{'time',480},'Format','netcdf4');
ncwrite('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','time',timenum);
ncwriteatt('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','time','standard_name','time');
ncwriteatt('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','time','long_name','time');
ncwriteatt('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','time','units','months since 1979-12-15 00:00:00');
ncwriteatt('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','time','calendar','standard');
ncwriteatt('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','time','axis','T');
ncwriteatt('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','time','cell_methods','time: mean');

nccreate('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','x','dimensions',{'x',720});
ncwrite('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','x',x);
ncwriteatt('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','x','standard_name','longitude');
ncwriteatt('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','x','long_name','longitude');
ncwriteatt('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','x','units','degrees_E');
ncwriteatt('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','x','axis','X');

nccreate('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','y','dimensions',{'y',360});
ncwrite('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','y',y);
ncwriteatt('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','y','standard_name','latitude');
ncwriteatt('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','y','long_name','latitude');
ncwriteatt('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','y','units','degrees_N');
ncwriteatt('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','y','axis','Y');

nccreate('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','discharge','Dimensions',{'time','y','x'},'ChunkSize',[480 360 1]);
ncwrite('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','discharge',gridded_discharg);
ncwriteatt('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','discharge','standard_name','discharge');
ncwriteatt('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','discharge','long_name','river_discharge');
ncwriteatt('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','discharge','units','m^3/s');
ncwriteatt('/aosc/horse/carton/rivers/gridded_discharge_mo_360x720x480R.nc','discharge','cell_methods','time: mean');

