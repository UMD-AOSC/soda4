function [Ba_total_discharge_80_19,Ba_lon,Ba_lat]= Add_Greenland
%
% Data from "Bamber, et al. (2018) Land ice freshwater budget of the Arctic and North Atlantic Oceans. 
% Part I: Data, methods and results. J. Geophys. Res.  https://doi.org/10.1002/2017JC013605"
% total should be: 1,300 km3/yr (2010-2016).  We sum ice solid discharge, runoff, + tundra discharges (giving 1336 km3/yr)
% X: 752, Y:785, time: 708
%ncdisp('/aosc/horse/carton/rivers/greenland/bamber18/FWF17.v3_b.nc')
%
%
Land = ncread('/aosc/horse/carton/rivers/greenland/bamber18/FWF17.v3_b.nc','LSMGr'); %Hole-filled Greenland land mass mask ocean=0, land=1
%contourf(flipud(Land'));colorbar;
Ba_lat = ncread('/aosc/horse/carton/rivers/greenland/bamber18/FWF17.v3_b.nc','lat'); % lat of Bamber discharge points [752x785]
Ba_lon = ncread('/aosc/horse/carton/rivers/greenland/bamber18/FWF17.v3_b.nc','lon'); % lon of Bamber discharge points [752x785]
jday = ncread('/homes/metofac/carton/transfer/FWF17.v3_b.nc','TIME')+juliandate(1958,01,01); % days since 1/1/1958 [708x1]
%
% 1) Ice sheet runoff just from Greenland [752x785x708] (the RGRIS term of Bamber (2018) eqn. 1)
% scale_factor = 0.01 units = 'km3' uncertainty: 20%
%
Ba_runoff_ice = ncread('/aosc/horse/carton/rivers/greenland/bamber18/FWF17.v3_b.nc','runoff_ice').*Land; 
%
% Check: plot the RGIS time series (2yr running average) to match Bamber et al. (2018) Fig. 3
%
Ba_runoff_ice_ts = squeeze(sum(Ba_runoff_ice,[1 2])); %units: km3/mo
%plot(movmean(Ba_runoff_ice_ts,[11 12])*12);%units: km3/yr 
%hold on;
%
% 2) Solid ice discharge just from Greenland [752x785x708] (the D term of Bamber (2018) eqn. 1) 
% scale_factor = 0.01 units = 'km3' uncertainty: 10%
%
Ba_solid_ice = ncread('/aosc/horse/carton/rivers/greenland/bamber18/FWF17.v3_b.nc','solid_ice').*Land;; 
%
% Check: plot the D time series (2yr running average) to match Bamber et al. (2018) Fig. 3
%
%Ba_solid_ice_ts = squeeze(sum(Ba_solid_ice,[1 2])); %units: km3/mo
%plot(movmean(Ba_solid_ice_ts,[11 12])*12); %units: km3/yr
%
% 3) Tundra runoff [752x785x708] (the Rt term of Bamber (2018) eqn. 1) 
% scale_factor = 0.01 units = 'km3' uncertainty: 10%
%
Ba_runoff_tundra = ncread('/aosc/horse/carton/rivers/greenland/bamber18/FWF17.v3_b.nc','runoff_tundra').*Land;; 
%
% Check: plot the Rt time series (2yr running average) to match Bamber et al. (2018) Fig. 3
%
Ba_runoff_tundra_ts = squeeze(sum(Ba_runoff_tundra,[1 2])); %units: km3/mo
%plot(movmean(Ba_runoff_tundra_ts,[11 12])*12); %units: km3/yr 
%hold off;
%
% combine the three contributions and convert to m3/s
%
Ba_total_discharge_80_17 = (Ba_runoff_ice(:,:,253:708) + Ba_solid_ice(:,:,253:708) + Ba_runoff_tundra(:,:,253:708))*380; %units: m3/s (mult by 380)
%
Ba_total_discharge_16 = Ba_total_discharge_80_17(:,:,433:444);
%
% extend the data set by adding two copies of year 2016 after year 2017 (I trust 2016 a bit more than 2017)
%
Ba_total_discharge_80_19 = cat(3,Ba_total_discharge_80_17,Ba_total_discharge_16,Ba_total_discharge_16);
%
% Combined monthly file: Ba_total_discharge_80_19      [752x785x480] units: m3/s
%plot(movmean(squeeze(sum(Ba_total_discharge_80_19,[1 2])),[2 3])); title('Total Greenland Discharge [m3/s]');
end
