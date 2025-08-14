%
% Greenland discharge estimates from 
% Mankoff et al., 2019: https://doi.org/10.5194/essd-11-769-2019
%
lat = ncread('gate.nc','mean_lat'); % mean lat of discharge gate [267x1]
lon = ncread('gate.nc','mean_lon'); % mean lon of discharge gate [267x1]
jday = ncread('gate.nc','time')+juliandate(1986,04,15); % days since 1986-04-15 [445x1]
discharge = ncread('gate.nc','discharge'); % land_ice_volume_flow_rate_due_to_calving_and_ice_front_melting [km3/yr] [445x267]
total_discharge = sum(discharge,2)*1.e9/(3.15e7*1.e6);
plot(total_discharge);
title('M20 TOTAL ICE DISCARGE');
ylabel('[Sv]');
xlabel('Months since 3/86');
datetime(jday(1:20),'convertfrom','juliandate')

