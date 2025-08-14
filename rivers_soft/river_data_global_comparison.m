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
%
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
f1 = figure('Name','Observations');  % global observations picture
f1.Position = [ 20 500 2000 800];
ts0 = timeseries(obs0);
ts0.TimeInfo.Format = 'mm/yyyy';
ts0.TimeInfo.Units = 'months';
ts0.TimeInfo.StartDate = '15-Jan-1900';
original_obs = plot(ts0);           % ds551.0 total river stations
original_obs.LineWidth = 1.5;
original_obs.Color = [ 0 1 0];
original_obs.DisplayName = 'original';
hold on
ts1 = timeseries(obs1);
ts1.TimeInfo.Format = 'mm/yyyy';
ts1.TimeInfo.Units = 'months';
ts1.TimeInfo.StartDate = '15-Jan-1900';
extend_obs = plot(ts1);             % original obs extended to 2023 and few rivers replaced by newer data
extend_obs.LineWidth = 1.5;
extend_obs.Color = [ 0 0 1];
extend_obs.DisplayName = 'extended';
ts2 = timeseries(obs2);
ts2.TimeInfo.Format = 'mm/yyyy';
ts2.TimeInfo.Units = 'months';
ts2.TimeInfo.StartDate = '15-Jan-1900';
filled_obs = plot(ts2);             % extended obs with gaps filled by climatology
filled_obs.LineWidth = 1.5;
filled_obs.Color = [ 1 0 0];
filled_obs.DisplayName = 'extended';
ax = gca;
ax.Title.String = 'Global river discharge observations';
ax.XLabel.String = 'Time';
ax.YLabel.String = 'Number of observations per month';
%
lg = legend('ds551.0','extended','filled');
lg.Position = [ 0.80 0.75 0.075 0.05];
hold off
%
%
f2 = figure('Name','Global river discharge');  % global river discharge picture
f2.Position = [ 100 400 2000 800];
tsg0 = timeseries(glob0c);                     % global river discharge from ds551.0
tsg0.TimeInfo.Format = 'mm/yyyy';
tsg0.TimeInfo.Units = 'months';
tsg0.TimeInfo.StartDate = '15-Jan-1900';
original_dis = plot(tsg0);     
original_dis.LineWidth = 1.5;                
original_dis.Color = [ 0 1 0];
original_dis.DisplayName = 'original';
hold on
tsg1 = timeseries(glob1c);                     % global river discharge from extended ds551.0
tsg1.TimeInfo.Format = 'mm/yyyy';
tsg1.TimeInfo.Units = 'months';
tsg1.TimeInfo.StartDate = '15-Jan-1900';
extend_dis = plot(tsg1);     
extend_dis.LineWidth = 1.5;                
extend_dis.Color = [ 0 0 1];
extend_dis.DisplayName = 'extended';
tsg2 = timeseries(glob2c);                     % global river discharge from extended ds551.0 with filled gaps
tsg2.TimeInfo.Format = 'mm/yyyy';
tsg2.TimeInfo.Units = 'months';
tsg2.TimeInfo.StartDate = '15-Jan-1900';
extend_dis = plot(tsg2);     
extend_dis.LineWidth = 1.5;                
extend_dis.Color = [ 1 0 0];
extend_dis.DisplayName = 'extended';
ax2 = gca;
ax2.Title.String = 'Global river discharge';
ax2.XLabel.String = 'Time';
ax2.YLabel.String = 'Global river discharge (m3/sec)';
%
mean0 = 'ds551.0 (mean '+string(mean(glob0c,'omitnan'))+' m3/sec)';
mean1 = 'extended (mean '+string(mean(glob1c,'omitnan'))+' m3/sec)';
mean2 = 'filled (mean '+string(mean(glob2c,'omitnan'))+' m3/sec)';
lg = legend(mean0,mean1,mean2);
lg.Position = [ 0.65 0.15 0.15 0.05];
hold off
%
%
f3 = figure('Name','Oceans river discharge');  % different oceans river discharge picture (from filled data)
f3.Position = [ 200 300 2000 800];
tsat = timeseries(atl2c+med2c);                % Atlantic ocean + Mediterranean sea river discharge
tsat.TimeInfo.Format = 'mm/yyyy';
tsat.TimeInfo.Units = 'months';
tsat.TimeInfo.StartDate = '15-Jan-1900';
atl_dis = plot(tsat);     
atl_dis.LineWidth = 1.5;                
atl_dis.Color = [ 0 0 1];
atl_dis.DisplayName = 'Atlantic';
hold on
tsar = timeseries(arc2c);                      % Arctic ocean river discharge
tsar.TimeInfo.Format = 'mm/yyyy';
tsar.TimeInfo.Units = 'months';
tsar.TimeInfo.StartDate = '15-Jan-1900';
arc_dis = plot(tsar);     
arc_dis.LineWidth = 1.5;                
arc_dis.Color = [ 1 0 0];
arc_dis.DisplayName = 'Arctic';
tsin = timeseries(ind2c);                      % Indian ocean river discharge
tsin.TimeInfo.Format = 'mm/yyyy';
tsin.TimeInfo.Units = 'months';
tsin.TimeInfo.StartDate = '15-Jan-1900';
ind_dis = plot(tsin);     
ind_dis.LineWidth = 1.5;                
ind_dis.Color = [ 0 1 0];
ind_dis.DisplayName = 'Indian';
tspa = timeseries(pac2c);                      % Pacific ocean river discharge
tspa.TimeInfo.Format = 'mm/yyyy';
tspa.TimeInfo.Units = 'months';
tspa.TimeInfo.StartDate = '15-Jan-1900';
pac_dis = plot(tspa);     
pac_dis.LineWidth = 1.5;                
pac_dis.Color = [ 1 1 0];
pac_dis.DisplayName = 'Pacific';
%
ax3 = gca;
ax3.Title.String = 'Oceans river discharge';
ax3.XLabel.String = 'Time';
ax3.YLabel.String = 'River discharge (m3/sec)';
%
meanat = 'Atlantic (mean '+string(mean(atl2c+med2c,'omitnan'))+' m3/sec)';
meanar = 'Arctic (mean '+string(mean(arc2c,'omitnan'))+' m3/sec)';
meanin = 'Indian (mean '+string(mean(ind2c,'omitnan'))+' m3/sec)';
meanpa = 'Pacific (mean '+string(mean(pac2c,'omitnan'))+' m3/sec)';
lg = legend(meanat,meanar,meanin,meanpa);
lg.Position = [ 0.15 0.85 0.15 0.05];
hold off
%
% Adjust rivers discharge to the unmonitored areas in each ocean follow to 
% Dai and Trenberth (2002) (exclude 10 stronges rivers)
flow3 = flow2;
nn = 3;
ind = indatl(find(indatl>nn)); flow3(ind,:) = flow2(ind,:)*1.2793;  % adjust atlantic ocean rivers
ind = indarc(find(indarc>nn)); flow3(ind,:) = flow2(ind,:)*1.3479;  % adjust arctic ocean rivers
ind = indind(find(indind>nn)); flow3(ind,:) = flow2(ind,:)*1.6946;  % adjust indian ocean rivers
ind = indpac(find(indpac>nn)); flow3(ind,:) = flow2(ind,:)*2.2305;  % adjust pacific ocean rivers
%
atl3 = sum(flow3(indatl,:).*ratio_m2s(indatl),1,'omitnan');  % Atlantic ocean rivers discharge
arc3 = sum(flow3(indarc,:).*ratio_m2s(indarc),1,'omitnan');  % Arctic ocean rivers discharge
ind3 = sum(flow3(indind,:).*ratio_m2s(indind),1,'omitnan');  % Indian ocean rivers discharge
med3 = sum(flow3(indmed,:).*ratio_m2s(indmed),1,'omitnan');  % Mediterranean sea rivers discharge
pac3 = sum(flow3(indpac,:).*ratio_m2s(indpac),1,'omitnan');  % Pacific ocean rivers discharge
%
meanat3 = 'Atlantic (mean '+string(mean(atl3+med3,'omitnan'))+' m3/sec)';
meanar3 = 'Arctic (mean '+string(mean(arc3,'omitnan'))+' m3/sec)';
meanin3 = 'Indian (mean '+string(mean(ind3,'omitnan'))+' m3/sec)';
meanpa3 = 'Pacific (mean '+string(mean(pac3,'omitnan'))+' m3/sec)';
%
f4 = figure('Name','Oceans river discharge adjusted');  % Adjusted for each ocean river discharge picture (from filled data)
f4.Position = [ 300 200 2000 800];
tsat3 = timeseries(atl3+med3);                % Atlantic ocean + Mediterranean sea river discharge
tsat3.TimeInfo.Format = 'mm/yyyy';
tsat3.TimeInfo.Units = 'months';
tsat3.TimeInfo.StartDate = '15-Jan-1900';
at3_dis = plot(tsat3);     
at3_dis.LineWidth = 1.5;                
at3_dis.Color = [ 0 0 1];
at3_dis.DisplayName = 'Atlantic';
hold on
tsar3 = timeseries(arc3);                     % Arctic ocean river discharge
tsar3.TimeInfo.Format = 'mm/yyyy';
tsar3.TimeInfo.Units = 'months';
tsar3.TimeInfo.StartDate = '15-Jan-1900';
ar3_dis = plot(tsar3);     
ar3_dis.LineWidth = 1.5;                
ar3_dis.Color = [ 1 0 0];
ar3_dis.DisplayName = 'Arctic';
tsin3 = timeseries(ind3);                      % Indian ocean river discharge
tsin3.TimeInfo.Format = 'mm/yyyy';
tsin3.TimeInfo.Units = 'months';
tsin3.TimeInfo.StartDate = '15-Jan-1900';
in3_dis = plot(tsin3);     
in3_dis.LineWidth = 1.5;                
in3_dis.Color = [ 0 1 0];
in3_dis.DisplayName = 'Indian';
tspa3 = timeseries(pac3);                      % Pacific ocean river discharge
tspa3.TimeInfo.Format = 'mm/yyyy';
tspa3.TimeInfo.Units = 'months';
tspa3.TimeInfo.StartDate = '15-Jan-1900';
pa3_dis = plot(tspa);     
pa3_dis.LineWidth = 1.5;                
pa3_dis.Color = [ 1 1 0];
pa3_dis.DisplayName = 'Pacific';
%
ax4 = gca;
ax4.Title.String = 'Oceans river discharge (adjusted for each ocean)';
ax4.XLabel.String = 'Time';
ax4.YLabel.String = 'River discharge (m3/sec)';
%
meanat3 = 'Atlantic (mean '+string(mean(atl3+med3,'omitnan'))+' m3/sec)';
meanar3 = 'Arctic (mean '+string(mean(arc3,'omitnan'))+' m3/sec)';
meanin3 = 'Indian (mean '+string(mean(ind3,'omitnan'))+' m3/sec)';
meanpa3 = 'Pacific (mean '+string(mean(pac3,'omitnan'))+' m3/sec)';
lg = legend(meanat3,meanar3,meanin3,meanpa3);
lg.Position = [ 0.15 0.85 0.15 0.05];
hold off
%
