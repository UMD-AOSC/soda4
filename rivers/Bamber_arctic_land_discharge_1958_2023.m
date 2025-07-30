%
close all hidden; 
clearvars;
%
rdir = "/aosc/iceland2/chepurin/RIVERS/DATA/";    % data directory
fname = rdir+"GREENLAND/bamber18/FWF17.v3_b.nc";  % Greenland and tundra runoff data file
%
irunoff = ncread(fname,'runoff_ice');     % read Ice sheet runoff data
srunoff = ncread(fname,'solid_ice');      % read Solid ice discharge data
trunoff = ncread(fname,'runoff_tundra');  % read Tundra runoff data
glon = ncread(fname,'lon');               % read longitude data
glat = ncread(fname,'lat');               % read latitude data
timn = ncread(fname,'TIME');              % read time data
%
time = (floor(([1:792]-0.5)/12)+1958)*100+[1:792]-floor(([1:792]-0.5)/12)*12; % time in format YYYYMM
%
% convert grid data to station data and units to m3/sec
lon = reshape(glon,[],1);  % station longitude
lat = reshape(glat,[],1);  % station latitude
dum = reshape(irunoff,[],708);
dumind = find(sum(dum,2));
ir = dum(dumind,:)*380.5;  % ice sheet runoff data [m3/sec]
lon_ir = lon(dumind);      % ice sheet runoff data lonritude
lat_ir = lat(dumind);      % ice sheet runoff data latitude
dum = reshape(srunoff,[],708);
dumind = find(sum(dum,2));
is = dum(dumind,:)*380.5;  % solid ice discharge data [m3/sec]
ind = find(is<0.0);        % filter negative ice discharge data
is(ind) = 0.0;             % filter negative ice discharge data
lon_is = lon(dumind);      %  solid ice discharge data lonritude
lat_is = lat(dumind);      %  solid ice discharge data latitude
dum = reshape(trunoff,[],708);
dumind = find(sum(dum,2));
it = dum(dumind,:)*380.5;  % tundra runoff data [m3/sec]          
lon_it = lon(dumind);      % tundra runoff data lonritude
lat_it = lat(dumind);      % tundra runoff data latitude
%clear dum dumind irunoff trunoff srunoff glon glat;
%
% fille gaps and extend data by climatology (2006-2016) [589:708]
irex(:,1:708) = ir;
itex(:,1:708) = it;
isex(:,1:708) = is;
for i = 1:size(ir,1);   % loop over stations
    irex(i,709:792) = repmat(mean(reshape(ir(i,589:708),12,[]),2,'omitmissing'),7,1);
end
for i = 1:size(it,1);   % loop over stations
    itex(i,709:792) = repmat(mean(reshape(it(i,589:708),12,[]),2,'omitmissing'),7,1);
end
for i = 1:size(is,1);   % loop over stations
    isex(i,709:792) = repmat(mean(reshape(is(i,589:708),12,[]),2,'omitmissing'),7,1);
end
%
tir = sum(irex,1);  % total Greenland ice sheet runoff
tis = sum(isex,1);  % total Greenland solid ice discharge
tit = sum(itex,1);  % total arctic tundra runoff
%
% draw Arctic land fresh water total fluxes picture
f1 = figure('Name','Arctic fresh water fluxes from the land');   % total Arctic land discharge picture
f1.Position = [ 20 500 2000 1300];
%
ts1 = timeseries(tir+tis+tit);
ts1.TimeInfo.Format = 'mm/yyyy';
ts1.TimeInfo.Units = 'months';
ts1.TimeInfo.StartDate = '15-Jan-1958';
mean1 = 'total discharge mean '+string(mean(tir+tis+tit,'omitnan'))+' m3/sec)';
pic1 = plot(ts1);
pic1.LineWidth = 2.0;
pic1.LineStyle = '-';
pic1.Color = [ 0 0 0];
pic1.DisplayName = 'ice sheet runoff';
hold on
%
ts2 = timeseries(tis);
ts2.TimeInfo.Format = 'mm/yyyy';
ts2.TimeInfo.Units = 'months';
ts2.TimeInfo.StartDate = '15-Jan-1958';
mean2 = 'Greenland solid ice discharge mean '+string(mean(tis,'omitnan'))+' m3/sec)';
pic2 = plot(ts2);
pic2.LineWidth = 2.0;
pic2.LineStyle = '-';
pic2.Color = [ 0 1 0];
%
ts3 = timeseries(tir);
ts3.TimeInfo.Format = 'mm/yyyy';
ts3.TimeInfo.Units = 'months';
ts3.TimeInfo.StartDate = '15-Jan-1958';
mean3 = 'Greenland ice sheet runoff mean '+string(mean(tir,'omitnan'))+' m3/sec)';
pic3 = plot(ts3);
pic3.LineWidth = 2.0;
pic3.LineStyle = '-';
pic3.Color = [ 1 0 0];
%
ts4 = timeseries(tit);
ts4.TimeInfo.Format = 'mm/yyyy';
ts4.TimeInfo.Units = 'months';
ts4.TimeInfo.StartDate = '15-Jan-1958';
mean4 = 'Tundra runoff mean '+string(mean(tit,'omitnan'))+' m3/sec)';
pic4 = plot(ts4);
pic4.LineWidth = 2.0;
pic4.LineStyle = '-';
pic4.Color = [ 0 0 1];
%
ax = gca;
ax.Title.String = 'Arctic fresh water fluxes from the land';
ax.TitleFontSizeMultiplier = 2.5;
ax.TitleFontWeight = 'normal';
ax.XLabel.String = 'Time';
ax.YLabel.String = 'Fresh water flux (m3/sec)';
%
lg = legend(mean1,mean2,mean3,mean4);
lg.Position = [ 0.16 0.85 0.15 0.06];
%

% write extended discharge data in NetCDF file
fnameout = rdir+'Bamber_arctic_land_discharge_ext_1958_2023.nc';
nccreate(fnameout,'time','dimension',{'time',Inf},'Datatype','int32','Format','netcdf4');
ncwrite(fnameout,'time',time);
ncwriteatt(fnameout,'time','long_name','time as YYYYMM');
ncwriteatt(fnameout,'time','units','unitless');
%
[dm,d] = size(irex);
st_ir = [0:dm-1];      %  
nccreate(fnameout,'lon_ir','dimension',{'st_ir',dm},'Datatype','single','Format','netcdf4');
ncwrite(fnameout,'lon_ir',lon_ir);
ncwriteatt(fnameout,'lon_ir','long_name','ice sheet runoff data longitude');
ncwriteatt(fnameout,'lon_ir','units','degrees_east');
nccreate(fnameout,'lat_ir','dimension',{'st_ir',dm},'Datatype','single','Format','netcdf4');
ncwrite(fnameout,'lat_ir',lat_ir);
ncwriteatt(fnameout,'lat_ir','long_name','ice sheet runoff data latitude');
ncwriteatt(fnameout,'lat_ir','units','degrees_north');
nccreate(fnameout,'ice_runoff','dimension',{'st_ir',dm,'time',792},'FillValue',NaN,'Datatype','single','Format','netcdf4');
ncwrite(fnameout,'ice_runoff',irex)  % ice runoff data extended to 12/2023
ncwriteatt(fnameout,'ice_runoff','long_name','Ice sheet runoff monthly mean values');
ncwriteatt(fnameout,'ice_runoff','units','m3/sec');
%
[dm,d] = size(itex);
st_it = [0:dm-1];      %  
nccreate(fnameout,'lon_it','dimension',{'st_it',dm},'Datatype','single','Format','netcdf4');
ncwrite(fnameout,'lon_it',lon_it);
ncwriteatt(fnameout,'lon_it','long_name','tundra runoff data longitude');
ncwriteatt(fnameout,'lon_it','units','degrees_east');
nccreate(fnameout,'lat_it','dimension',{'st_it',dm},'Datatype','single','Format','netcdf4');
ncwrite(fnameout,'lat_it',lat_it);
ncwriteatt(fnameout,'lat_it','long_name','tundra runoff data latitude');
ncwriteatt(fnameout,'lat_it','units','degrees_north');
nccreate(fnameout,'tundra_runoff','dimension',{'st_it',dm,'time',792},'FillValue',NaN,'Datatype','single','Format','netcdf4');
ncwrite(fnameout,'tundra_runoff',itex)  % tundra runoff data extended to 12/2023
ncwriteatt(fnameout,'tundra_runoff','long_name','Tundra runoff monthly mean values');
ncwriteatt(fnameout,'tundra_runoff','units','m3/sec');
%
[dm,d] = size(isex);
st_is = [0:dm-1];      %  
nccreate(fnameout,'lon_is','dimension',{'st_is',dm},'Datatype','single','Format','netcdf4');
ncwrite(fnameout,'lon_is',lon_is);
ncwriteatt(fnameout,'lon_is','long_name','solid ice discharge data longitude');
ncwriteatt(fnameout,'lon_is','units','degrees_east');
nccreate(fnameout,'lat_is','dimension',{'st_is',dm},'Datatype','single','Format','netcdf4');
ncwrite(fnameout,'lat_is',lat_is);
ncwriteatt(fnameout,'lat_is','long_name','solid_ice discharge data latitude');
ncwriteatt(fnameout,'lat_is','units','degrees_north');
nccreate(fnameout,'solid_ice','dimension',{'st_is',dm,'time',792},'FillValue',NaN,'Datatype','single','Format','netcdf4');
ncwrite(fnameout,'solid_ice',isex)  % solid ice discharge data extended to 12/2023
ncwriteatt(fnameout,'solid_ice','long_name','Solid ice discharge monthly mean values');
ncwriteatt(fnameout,'solid_ice','units','m3/sec');
%



