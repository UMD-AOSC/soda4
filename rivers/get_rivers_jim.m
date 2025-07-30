function [fa_full,fa_clim,fa_lon,fa_lat] = get_rivers
%
% Load Dai-Trenberth river discharge data.  To check file content: ncdisp('/aosc/horse/carton/rivers/dai-tren/coastal-stns-Vol-monthly.updated-May2019.nc');
% units: m3/s
%
fa_lon = ncread('/aosc/horse/carton/rivers/dai-tren/coastal-stns-Vol-monthly.updated-May2019.nc','lon_mou'); % lon of river mouth [925x1]
fa_lat = ncread('/aosc/horse/carton/rivers/dai-tren/coastal-stns-Vol-monthly.updated-May2019.nc','lat_mou'); % lat of river mouth [925x1]
flow = ncread('/aosc/horse/carton/rivers/dai-tren/coastal-stns-Vol-monthly.updated-May2019.nc','FLOW'); % monthly mean transport for each river gauge. Begins 1/1/1900 [925x1428] [m3/s]
ratio_m2s = ncread('/aosc/horse/carton/rivers/dai-tren/coastal-stns-Vol-monthly.updated-May2019.nc','ratio_m2s'); % ratio of transport at river mouth and at gauge station [925x1]
riv_name = ncread('/aosc/horse/carton/rivers/dai-tren/coastal-stns-Vol-monthly.updated-May2019.nc','riv_name'); % river name [30x925]
ocn_name = ncread('/aosc/horse/carton/rivers/dai-tren/coastal-stns-Vol-monthly.updated-May2019.nc','ocn_name'); % name of ocean into which river discharges
time = ncread('/aosc/horse/carton/rivers/dai-tren/coastal-stns-Vol-monthly.updated-May2019.nc','time'); % 'time as YYYYMM' [1428x1] [961 is 1/1/1980]
%
% update delta location for 15 key rivers (largest/Arctic)
%
fa_lat(1,1)=0.455;%Amazon
fa_lon(1,1)=-49.86;
fa_lat(2,1)=-6.04;%Congo
fa_lon(2,1)=12.34;
fa_lat(3,1)=9.00;%Orinoco
fa_lon(3,1)=-60.72;
fa_lat(4,1)=31.33;%Yangtze
fa_lon(4,1)=121.94;
fa_lat(5,1)=22.60;%Ganges
fa_lon(5,1)=90.88;
fa_lat(6,1)=29.15;%Mississippi
fa_lon(6,1)=-89.24;
fa_lat(7,1)=72.56;%Lena
fa_lon(7,1)=79.88;
fa_lat(8,1)=-34.30;%Plata
fa_lon(8,1)=-58.35;
fa_lat(9,1)=72.87;%Lena
fa_lon(9,1)=129.73;
fa_lat(13,1)=66.47;%Ob
fa_lon(13,1)=71.52;
fa_lat(19,1)=69.58;%Mackenzie
fa_lon(19,1)=-133.66;
fa_lat(24,1)=62.58;%Yukon
fa_lon(24,1)=-164.92;
fa_lat(31,1)=68.23;%Pechora
fa_lon(31,1)=54.25;
fa_lat(33,1)=74.14;%Khatanga
fa_lon(33,1)=110.36;
fa_lat(35,1)=69.59;%Kolyma
fa_lon(35,1)=161.31;

%
% Compute the monthly climatology for each river.  We'll use this to fill in missing values
%
fa_trans = flow(:,961:1428); %load the 40-yr monthly discharge 1/1980-12/2018
fa1 = fa_trans';%put flow into matrix fa with time in rows, rivers in columns [468,925]
for irev = 1:size(fa1,2); % correct for the fact that the measurements are made at gauge, not river mouth
	fa(:,irev)=fa1(:,irev)*ratio_m2s(irev);
end
%
% Check global average transport: mean(nansum(fa,2)) = 6.1665e+05 m3/s
%The total river discharge as reported by Dai and Trenberth (2002) is 1.18Sv.  
clear fa_trans fa1
%
fa_clim = zeros(12,size(fa,2)); %initialize the climatology [12x925]
for irev = 1:size(fa,2); % loop over 925 rivers and reshape each so that there are 12 columns for the 12  months
	far = reshape(fa(:,irev),12,[]);
	faa = nanmean(far,2); %compute the climotological monthly average ignoring NANs
	fa_clim(:,irev) = faa; %load it into fa_clim
end
clear far faa
%
% Check an example river: fa_clim(:,1) is the Amazon.  Annual avg according to Wikipedia: 2.09x105 m3/s.  Here we get 2.16x105 m3/s
%
% build the 40-yr long river discharge array fa_full by beginning with monthly climatology and inserting observed values
%
fa_full = repmat(fa_clim,40,1); % initialize fa_full with the monthly climatology 1/1980-12/2019 [480,925]
for irev = 1:size(fa,2)
	for itime = 1:468
		if ~isnan(fa(itime,irev)) 
			fa_full(itime,irev)=fa(itime,irev); % insert observed values
end;end;end;
%
%  Check Orinoco river transport to see the impact.  hold off;plot(fa_full(:,3));hold on;plot(fa(:,3));title('ORINOCO TRANSPORT');
%  Check Lena river transport to see the impact.  hold off;plot(fa_full(:,9));hold on;plot(fa(:,9));title('LENA TRANSPORT');
%  Check Ob river transport to see the impact.  hold off;plot(fa_full(:,13));hold on;plot(fa(:,13));title('Ob TRANSPORT');

%
%  Check mean global river discharge: Mean(nansum(fa_full,2))  = 0.803Sv is 32% less than the value of 1.18Sv given in Dai and Trenberth (2002)
%  Note: most of this difference is due to basin-dependent inflation of transports of ~30% by Dai and Trenberth (2002).  I carry out this inflation in 
%  a separate GRADS .exec (gridded_discharge_mo_360x720x480.exec)
end 
