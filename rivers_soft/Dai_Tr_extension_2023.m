%
close all hidden
%
% river_fname is a string vector which contains file's names with new data
% of rivers to be modified in the original Dai_Trenberth data file
% (coastal-stns-Vol-monthly.updated-May2019.nc)
rdir = "/aosc/iceland2/chepurin/RIVERS/DATA/";  % dirrectory where are the river discharge data
river_fname(1) = rdir+"Amazon_Obidos_discharge_1968_2020.txt";
river_fname(2) = rdir+"Congo_kinshasa_discharge_1903-2010.txt";
river_fname(3) = rdir+"Orinoco_Ptuente_Angostura_discharge_1923-1989.txt";
river_fname(4) = rdir+"Mississippi_Vicksburg_discharge_1928_2022.txt";
river_fname(5) = rdir+"Parana_Timmbues_discharge_1905_2021.txt";
river_fname(6) = rdir+"Irrawaddy_Sagaing_discharge_1978_2015.txt";
river_fname(7) = rdir+"St_Lawrence_Cornwall_discharge_1935-2023.txt";
river_fname(8) = rdir+"Xingu_Altamira_discharge_1968_2020.txt";
river_fname(9) = rdir+"Yenisei_Igarka_mn_1936-2023.txt";
river_fname(10) = rdir+"Lena_Kyusyur_mn_1936-2023.txt";
river_fname(11) = rdir+"Ob_Salekhard_mn_1936-2023.txt";
river_fname(12) = rdir+"Mackenzie_ArcticRedRiver_mn_1972-2023.txt";
river_fname(13) = rdir+"Yukon_PilotStation_mn_1975-2023.txt";
river_fname(14) = rdir+"Pechora_UstTsilma_mn_1932-2023.txt";
river_fname(15) = rdir+"Kolyma_Kolymskoe_mn_1978-2023.txt";
river_fname(16) = rdir+"NorthernDvina_UstPenega_mn_1936-2023.txt";
river_fname(17) = rdir+"Olenek_7.5_km_down_of_Buur.s_mouth_mn_1991-2023.txt";
river_fname(18) = rdir+"Mezen_Malonisogorskoe_mn_1978-2023.txt";
river_fname(19) = rdir+"Nadym_Nadym_mn_1978-2023.txt";
river_fname(20) = rdir+"Onega_Porog_mn_1978-2023.txt";
%
% river_nt: 
% 1 column - river order in the original Dai_Trenberth's data file
% 2 column - month's number relative to Jan 1900 data modificateon starts
% 3 column - month's number relative to Jan 1900 data modificateon stopted
river_nt(1,:) = [1,817,1441];     % Amazon
river_nt(2,:) = [2,37,1332];      % Congo
river_nt(3,:) = [3,281,1080];     % Orinoco
river_nt(4,:) = [6,337,1473];     % Mississippi
river_nt(5,:) = [8,61,1460];      % Parana
river_nt(6,:) = [15,937,1392];    % Irrawaddy
river_nt(7,:) = [16,430,1477];    % St Lawrence
river_nt(8,:) = [18,822,1446];    % Xingu (Brazil, low Amazon)
river_nt(9,:) = [7,433,1483];     % Yenisey
river_nt(10,:) = [9,433,1483];    % Lena
river_nt(11,:) = [13,433,1483];   % Ob
river_nt(12,:) = [19,867,1483];   % Mackenzie
river_nt(13,:) = [24,901,1479];   % Yukon
river_nt(14,:) = [31,385,1483];   % Pechora
river_nt(15,:) = [35,937,1483];   % Kolyma
river_nt(16,:) = [37,433,1483];   % Severnaya Dvina
river_nt(17,:) = [97,1093,1483];  % Olenek
river_nt(18,:) = [228,937,1483];  % Mezen
river_nt(19,:) = [239,937,1483];  % Nadym
river_nt(20,:) = [255,937,1483];  % Onega
%
% read river discharge data from separate files to the new_data
new_data(1:20,1:1488) = NaN;
for i = 1:8;   % read data from CRDC files
    new_data(i,river_nt(i,2):river_nt(i,3)) = read_crdc_river(river_fname(i));
end;
%
for i =9:20;   % read data from Arctic files
    new_data(i,river_nt(i,2):river_nt(i,3)) = read_Arctic_river(river_fname(i));
end;
%
get_river;    % read original Dai Trenberth data (NCAR Research Data Archive ds551.0)
%
flow(:,1429:1488) = NaN;   % extend original array to 07-2023
%
flow_new = flow;
for i =1:20   
   ind = find(isnan(new_data(i,:)) == 0);      % find cell indexes where new data exist
   flow_new(river_nt(i,1),ind) = new_data(i,ind);  % replace and extend original data 
end

% write extended discharge data in NetCDF file
fname = rdir+'Dai_Tr_extension_2023.nc';
time(1429:1488)=(floor(([1429:1488]-0.5)/12)+1900)*100 +[1429:1488]-floor(([1429:1488]-0.5)/12)*12;  % extend time to 06-2023
nccreate(fname,'time','dimension',{'time',Inf},'Datatype','int32','Format','netcdf4');
ncwrite(fname,'time',time);
ncwriteatt(fname,'time','long_name','time as YYYYMM');
ncwriteatt(fname,'time','units','unitless');
%
[dm,d] = size(stn);
station = [0:dm-1];
nccreate(fname,'station','dimension',{'station',dm},'Datatype','int32','Format','netcdf4');
ncwrite(fname,'station',station);
ncwriteatt(fname,'station','long_name','station index');
ncwriteatt(fname,'station','units','unitless');
%
nccreate(fname,'lon','dimension',{'station',dm},'Datatype','single','Format','netcdf4');
ncwrite(fname,'lon',lon);
ncwriteatt(fname,'lon','long_name','station longitude');
ncwriteatt(fname,'lon','units','degrees_east');
nccreate(fname,'lat','dimension',{'station',dm},'Datatype','single','Format','netcdf4');
ncwrite(fname,'lat',lat);
ncwriteatt(fname,'lat','long_name','station latitude');
ncwriteatt(fname,'lat','units','degrees_north');
nccreate(fname,'lon_mou','dimension',{'station',dm},'Datatype','single','Format','netcdf4');
ncwrite(fname,'lon_mou',lon_mou);
ncwriteatt(fname,'lon_mou','long_name','river mouth  longitude');
ncwriteatt(fname,'lon_mou','units','degrees_east');
nccreate(fname,'lat_mou','dimension',{'station',dm},'Datatype','single','Format','netcdf4');
ncwrite(fname,'lat_mou',lat_mou);
ncwriteatt(fname,'lat_mou','long_name','river mouth  latitude');
ncwriteatt(fname,'lat_mou','units','degrees_north');
%
nccreate(fname,'area_stn','dimension',{'station',dm},'Datatype','single','Format','netcdf4');
ncwrite(fname,'area_stn',area_stn');
ncwriteatt(fname,'area_stn','long_name','drainage areas at station');
ncwriteatt(fname,'area_stn','units','km2');
nccreate(fname,'area_mou','dimension',{'station',dm},'Datatype','single','Format','netcdf4');
ncwrite(fname,'area_mou',area_mou');
ncwriteatt(fname,'area_mou','long_name','drainage areas at river mouth ');
ncwriteatt(fname,'area_mou','units','km2');
nccreate(fname,'vol_stn','dimension',{'station',dm},'Datatype','single','Format','netcdf4');
ncwrite(fname,'vol_stn',vol_stn);
ncwriteatt(fname,'vol_stn','long_name','annual mean volume at station');
ncwriteatt(fname,'vol_stn','units','km3/yr');
nccreate(fname,'ratio_m2s','dimension',{'station',dm},'Datatype','single','Format','netcdf4');
ncwrite(fname,'ratio_m2s',ratio_m2s);
ncwriteatt(fname,'ratio_m2s','long_name','ratio of volume between river mouth and station');
ncwriteatt(fname,'ratio_m2s','units','unitless');
nccreate(fname,'elev','dimension',{'station',dm},'Datatype','single','Format','netcdf4');
ncwrite(fname,'elev',elev);
ncwriteatt(fname,'elev','long_name','elevation at station');
ncwriteatt(fname,'elev','units','m');
%
nccreate(fname,'ct_name','dimension',{'char',d,'station',dm},'Datatype','char','Format','netcdf4');
ncwrite(fname,'ct_name',ct_name);
ncwriteatt(fname,'ct_name','long_name','country of station');
nccreate(fname,'cn_name','dimension',{'char',d,'station',dm},'Datatype','char','Format','netcdf4');
ncwrite(fname,'cn_name',cn_name);
ncwriteatt(fname,'cn_name','long_name','continent of station');
nccreate(fname,'ocn_name','dimension',{'char',d,'station',dm},'Datatype','char','Format','netcdf4');
ncwrite(fname,'ocn_name',ocn_name);
ncwriteatt(fname,'ocn_name','long_name','ocean name of river discharge');
%
nccreate(fname,'stn_name','dimension',{'char',d,'station',dm},'Datatype','char','Format','netcdf4');
ncwrite(fname,'stn_name',stn');
ncwriteatt(fname,'stn_name','long_name','station name');
nccreate(fname,'riv_name','dimension',{'char',d,'station',dm},'Datatype','char','Format','netcdf4');
ncwrite(fname,'riv_name',river');
ncwriteatt(fname,'riv_name','long_name','river name');
%
nccreate(fname,'FLOW','dimension',{'station',dm,'time',1488},'FillValue',NaN,'Datatype','single','Format','netcdf4');
ncwrite(fname,'FLOW',flow_new);
ncwriteatt(fname,'FLOW','long_name','monthly mean volume at station');
ncwriteatt(fname,'FLOW','units','m3/s');


%
