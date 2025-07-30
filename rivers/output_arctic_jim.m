% find out what's in the file: ncdisp('coastal-stns-Vol-monthly.updated-May2019.nc');
lon = ncread('coastal-stns-Vol-monthly.updated-May2019.nc','lon_mou'); % lon of river mouth [925x1]
lat = ncread('coastal-stns-Vol-monthly.updated-May2019.nc','lat_mou'); % lat of river mouth [925x1]
flow = ncread('coastal-stns-Vol-monthly.updated-May2019.nc','FLOW'); % monthly mean volume at station [925x1428] begins 1/1/1900 [1177 is 1/1/1980]
ratio_m2s = ncread('coastal-stns-Vol-monthly.updated-May2019.nc','ratio_m2s'); % ratio of volume between river mouth and station [925x1]
riv_name = ncread('coastal-stns-Vol-monthly.updated-May2019.nc','riv_name'); % river name [30x925]
ocn_name = ncread('coastal-stns-Vol-monthly.updated-May2019.nc','ocn_name'); % name of ocean into which river discharges
%time_1800 = ncread('coastal-stns-Vol-monthly.updated-May2019.nc','time_1800'); %  'months since 1800 01' [2412x1] [2161 is 1/1/1980] [ends 12/1/2000]
time = ncread('coastal-stns-Vol-monthly.updated-May2019.nc','time'); % 'time as YYYYMM' [1428x1] [961 is 1/1/1980]
7 Yenise 
9 Lena
13 Ob
19 Mackenzie
24 Yukon

fileID = fopen('yenise_1-1-1980.txt','w');
fprintf(fileID,'%8.0f\n',flow(7,1177:1428)*ratio_m2s(7));
fclose(fileID);
fileID = fopen('lena_1-1-1980.txt','w');
fprintf(fileID,'%8.0f\n',flow(9,1177:1428)*ratio_m2s(9));
fclose(fileID);
fileID = fopen('ob_1-1-1980.txt','w');
fprintf(fileID,'%8.0f\n',flow(13,1177:1428)*ratio_m2s(13));
fclose(fileID);
fileID = fopen('mackenzie_1-1-1980.txt','w');
fprintf(fileID,'%8.0f\n',flow(19,1177:1428)*ratio_m2s(19));
fclose(fileID);
fileID = fopen('yukon_1-1-1980.txt','w');
fprintf(fileID,'%8.0f\n',flow(24,1177:1428)*ratio_m2s(19));
fclose(fileID);
fileID = fopen('pechor_1-1-1980.txt','w');
fprintf(fileID,'%8.0f\n',flow(31,1177:1428)*ratio_m2s(19));
fclose(fileID);


