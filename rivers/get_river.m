close all hidden; clear everything

fname = '/aosc/horse/carton/rivers/dai-tren/coastal-stns-Vol-monthly.updated-May2019.nc';
river = (ncread(fname,'riv_name'))'; % river name [925x30]
stn = (ncread(fname,'stn_name'))'; % station name [925x30]
%stn_name = string(stn(:,1:30));
%for j = 1:925;
%    riv_name(j) = string(river(j,1));
%    for i = 2:30;
%        if river(j,i) == ' ';
%          break
%        end;
%    riv_name(j) = riv_name(j) + string(river(j,i));
%    end;
%end

ratio_m2s = ncread(fname,'ratio_m2s'); % ratio of volume between river mouth and station [925x1]
lon = ncread(fname,'lon'); % lon of station [925x1]
lat = ncread(fname,'lat'); % lat of station [925x1]
lon_mou = ncread(fname,'lon_mou'); % lon of river mouth [925x1]
lat_mou = ncread(fname,'lat_mou'); % lat of river mouth [925x1]
area_stn = ncread(fname,'area_stn'); % drainage areas at station [925x1]
area_mou = ncread(fname,'area_mou'); % drainage areas at river mouth [925x1]
vol_stn = ncread(fname,'vol_stn'); % annual volume at station [925x1]
elev = ncread(fname,'elev'); % elevation at station [925x1]
ct_name = ncread(fname,'ct_name'); % country of station [925x1]
cn_name = ncread(fname,'cn_name'); % continent of station [925x1]
ocn_name = ncread(fname,'ocn_name'); % ocean name of river discharge [925x1]
flow = ncread(fname,'FLOW'); % monthly mean volume at station [925x1428] begins 1/1/1900 [961 is 1/1/1980]
time = ncread(fname,'time'); % 'time as YYYYMM' [1428x1] [961 is 1/1/1980]

%k = 1;
%fileID = fopen(riv_name(k)+".txt",'w');
%fprintf(fileID,'%7s:',riv_name(k));
%formatSpec = ' lon %6.2f, lat %6.2f\n ratio = %6.4f \n';
%fprintf(fileID,formatSpec,lon(k),lat(k),ratio_m2s(k));
%fprintf(fileID,'%10.2f\n',flow(k,961:1428));
%fclose(fileID);
