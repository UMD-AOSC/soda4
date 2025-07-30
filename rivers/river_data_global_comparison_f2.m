%
close all hidden; 
clearvars;
%
rdir = "/aosc/iceland2/chepurin/RIVERS/DATA/";   % dirrectory where are the river discharge data
%
%fname0 = '/aosc/horse/carton/rivers/dai-tren/coastal-stns-Vol-monthly.updated-May2019.nc';
fname0 = rdir+'Dai_Tr_orig/coastal-stns-Vol-monthly.updated-May2019.nc'; % original Dai_Trenberth data (data set ds551.0 from NCAR archive)
fname1 = rdir+'Dai_Tr_extension_2023.nc';                    % fname0 updated by 20 rivers from GRDC and extendet to 2023
fname2 = rdir+'Dai_Tr_extension_filled_2023.nc';             % fname2 with gaps filled by climatology
%fname3 = rdir+'Dai_Tr_extension_filled_adjusted_ 2023.nc';  % fname2 with gaps filled by climatology and adgjusted for separate oceans
%
river = (ncread(fname0,'riv_name'))'; % river name [925x30]
stn = (ncread(fname0,'stn_name'))'; % station name [925x30]
ratio_m2s = ncread(fname0,'ratio_m2s'); % ratio of volume between river mouth and station [925x1]
lon = ncread(fname0,'lon'); % lon of station [925x1]
lat = ncread(fname0,'lat'); % lat of station [925x1]
lon_mou = ncread(fname0,'lon_mou'); % lon of river mouth [925x1]
lat_mou = ncread(fname0,'lat_mou'); % lat of river mouth [925x1]
area_stn = ncread(fname0,'area_stn'); % drainage areas at station [925x1]
area_mou = ncread(fname0,'area_mou'); % drainage areas at river mouth [925x1]
vol_stn = ncread(fname0,'vol_stn'); % annual volume at station [925x1]
elev = ncread(fname0,'elev'); % elevation at station [925x1]
ct_name = ncread(fname0,'ct_name'); % country of station [925x1]
cn_name = ncread(fname0,'cn_name'); % continent of station [925x1]
ocn_name = ncread(fname0,'ocn_name'); % ocean name of river discharge [925x1]
time = ncread(fname2,'time'); % 'time as YYYYMM' [1428x1] [961 is 1/1/1980]
%%
yyyy = floor(time/100);
mm = time-yyyy*100;
%
flow0 = ncread(fname0,'FLOW'); % monthly mean volume at station [925x1488] begins 1/1/1900 [961 is 1/1/1980]
flow1 = ncread(fname1,'FLOW'); % monthly mean volume at station [925x1488] begins 1/1/1900 [961 is 1/1/1980]
flow2 = ncread(fname2,'FLOW'); % monthly mean volume at station [925x1488] begins 1/1/1900 [961 is 1/1/1980]
%flow3 = ncread(fname3,'FLOW'); % monthly mean volume at station [925x1488] begins 1/1/1900 [961 is 1/1/1980]
%
obs0 = sum(~isnan(flow0(:,:)),1);  % total number of discharge observations per month from fname0
obs1 = sum(~isnan(flow1(:,:)),1);  % total number of discharge observations per month from fname1
obs2 = sum(~isnan(flow2(:,:)),1);  % total number of discharge observations per month from fname2
%
glob0 = sum(flow0,1,'omitnan');    % global monthly river discharge from fname0
glob1 = sum(flow1,1,'omitnan');    % global monthly river discharge from fname1
glob2 = sum(flow2,1,'omitnan');    % global monthly river discharge from fname2
%
glob0c = sum(flow0.*ratio_m2s,1,'omitnan'); % global monthly river discharge corrected for river station position from fname0
glob1c = sum(flow1.*ratio_m2s,1,'omitnan'); % global monthly river discharge corrected for river station position from fname1
glob2c = sum(flow2.*ratio_m2s,1,'omitnan'); % global monthly river discharge corrected for river station position from fname2
%
ocn = string(ocn_name(1:3,:)');    % get water basin name 
indatl = find(ocn == 'ATL');       % Atlantic ocean rivers
indpac = find(ocn == 'PAC');       % Pacific ocean rivers
indarc = find(ocn == 'ARC');       % Arctic ocean rivers
indind = find(ocn == 'IND');       % Indian ocean rivers
indmed = find(ocn == 'MED');       % Mediterranean sea rivers
%
atl2c = sum(flow2(indatl,:).*ratio_m2s(indatl),1,'omitnan');  % Atlantic ocean rivers discharge
arc2c = sum(flow2(indarc,:).*ratio_m2s(indarc),1,'omitnan');  % Arctic ocean rivers discharge
ind2c = sum(flow2(indind,:).*ratio_m2s(indind),1,'omitnan');  % Indian ocean rivers discharge
med2c = sum(flow2(indmed,:).*ratio_m2s(indmed),1,'omitnan');  % Mediterranean sea rivers discharge
pac2c = sum(flow2(indpac,:).*ratio_m2s(indpac),1,'omitnan');  % Pacific ocean rivers discharge
%
%
f2 = figure('Name','Global river discharge');  % global river discharge picture
f2.Position = [ 100 1500 2000 800];
tsg0 = timeseries(glob0c(960:1428));                     % global river discharge from ds551.0
tsg0.TimeInfo.Format = 'mm/yyyy';
tsg0.TimeInfo.Units = 'months';
tsg0.TimeInfo.StartDate = '15-Jan-1980';
original_dis = plot(tsg0);     
original_dis.LineWidth = 2.25;                
original_dis.Color = [ 1 0 0];
original_dis.DisplayName = 'original';
hold on
%tsg1 = timeseries(glob1c);                     % global river discharge from extended ds551.0
%tsg1.TimeInfo.Format = 'mm/yyyy';
%tsg1.TimeInfo.Units = 'months';
%tsg1.TimeInfo.StartDate = '15-Jan-1900';
%extend_dis = plot(tsg1);     
%extend_dis.LineWidth = 1.5;                
%extend_dis.Color = [ 0 1 0];
%extend_dis.DisplayName = 'extended';
tsg2 = timeseries(glob2c(960:1488));                     % global river discharge from extended ds551.0 with filled gaps
tsg2.TimeInfo.Format = 'mm/yyyy';
tsg2.TimeInfo.Units = 'months';
tsg2.TimeInfo.StartDate = '15-Jan-1980';
extend_dis = plot(tsg2);     
extend_dis.LineWidth = 2.25;                
extend_dis.Color = [ 0 0 1];
extend_dis.DisplayName = 'extended';
ax2 = gca;
ax2.FontSize = 16;
ax2.Title.String = 'Global discharge';
ax2.TitleFontWeight = 'normal';
ax2.Title.FontSize = 24;
ax2.XLabel.String = 'Time';
ax2.XLabel.FontSize = 18;
ax2.YLabel.String = 'Discharge (m3/sec)';
ax2.YLabel.FontSize = 16;
%
%mean0 = 'ds551.0 (mean '+string(mean(glob0c,'omitnan'))+' m3/sec)';
mean0 = 'Dai 2017 (mean = '+string(sprintf('%0.2e',mean(glob0c,'omitnan')))+' m3/sec)';
%mean1 = 'extended (mean '+string(mean(glob1c,'omitnan'))+' m3/sec)';
%mean2 = 'UMD (mean '+string(mean(glob2c,'omitnan'))+' m3/sec)';
mean2 = 'UMD (mean = '+string(sprintf('%0.2e',mean(glob2c,'omitnan')))+' m3/sec)';
%lg = legend(mean0,mean1,mean2);
lg = legend(mean0,mean2);
lg.FontSize = 18;
lg.Position = [ 0.18 0.15 0.15 0.05];
hold off
%
%