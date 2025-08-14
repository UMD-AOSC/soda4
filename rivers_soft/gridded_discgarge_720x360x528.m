%
close all hidden; 
clearvars;
%
rdir = "/aosc/iceland2/chepurin/RIVERS/DATA/";   % dirrectory where are the river discharge data
%
gfname = rdir+'Bamber_arctic_land_discharge_ext_1958_2023.nc';  % arctic land discharge data file
rfname = rdir+'Dai_Tr_extension_filled_adjusted_2023.nc';       % river discharge data file
%
discharge = zeros(720,360,528);   % initialize a 0.5x0.5 degree [90S-90N, 0E-360E] and monthly [1/1980-12/2023] array 
%
lon_ir = ncread(gfname,'lon_ir');  % read ice sheet runoff longitude
lat_ir = ncread(gfname,'lat_ir');  % read ice sheet runoff latitude
ir = ncread(gfname,'ice_runoff');  % read ice sheet runoff data
[ns,~] = size(ir);
for j =1:ns                    % loop over ice sheet runoff stations
  ix = round(lon_ir(j)*2)+361;  % convert station longitude to xgrid cell number
  iy = round(lat_ir(j)*2)+181;  % convert station latitude to ygrid cell number
  for mn = 1:528               % loop over months [starts from 1/1980 -> 264]
      discharge(ix,iy,mn) = discharge(ix,iy,mn)+ir(j,mn+264); 
  end
end
%
lon_it = ncread(gfname,'lon_it');     % read tundra runoff longitude
lat_it = ncread(gfname,'lat_it');     % read tundra runoff latitude
it = ncread(gfname,'tundra_runoff');  % read tundra runoff data
[ns,~] = size(it);
for j =1:ns                    % loop over ice sheet runoff stations
  ix = round(lon_it(j)*2)+361;  % convert station longitude to xgrid cell number
  iy = round(lat_it(j)*2)+181;  % convert station latitude to ygrid cell number
  for mn = 1:528               % loop over months [starts from 1/1980 -> 264]
      discharge(ix,iy,mn) = discharge(ix,iy,mn)+it(j,mn+264); 
  end
end
%
lon_is = ncread(gfname,'lon_is');     % read solid ice longitude
lat_is = ncread(gfname,'lat_is');     % read solid ice latitude
is = ncread(gfname,'solid_ice');      % read solid ice data
[ns,~] = size(is);
for j =1:ns                    % loop over solid ice stations
  ix = round(lon_is(j)*2)+361;  % convert station longitude to xgrid cell number
  iy = round(lat_is(j)*2)+181;  % convert station latitude to ygrid cell number
  for mn = 1:528               % loop over months [starts from 1/1980 -> 264]
      discharge(ix,iy,mn) = discharge(ix,iy,mn)+is(j,mn+264); 
  end
end
%
lon_dai = ncread(rfname,'lon_mou');     %  lon of river mouth [925x1]
lat_dai = ncread(rfname,'lat_mou');     %  lat of river mouth [925x1]
ratio_m2s = ncread(rfname,'ratio_m2s'); % ratio of volume between river mouth and station [925x1]
flow = ncread(rfname,'FLOW');  % monthly mean volume at station [925x1488] begins 1/1/1900 [961 is 1/1/1980]
time = ncread(rfname,'time');  % 'time as YYYYMM' [1488x1] [961 is 1/1/1980]
%
[ns,nt] = size(flow);
for j =1:ns                     % loop over ice sheet runoff stations
  ix = round(lon_dai(j)*2)+361;  % convert station longitude to xgrid cell number
  iy = round(lat_dai(j)*2)+181;  % convert station latitude to ygrid cell number
  for mn = 1:528                % loop over months [starts from 1/1980 -> 264]
      discharge(ix,iy,mn) = discharge(ix,iy,mn)+flow(j,mn+960)*ratio_m2s(j); 
  end
end
%
% write gridded discharge data in NetCDF file
tmn = 1:528;  % time starts at 01/1980 [month]
xlon = ((1:720)-361)/2; % grid cells longitude [180W-179.5E]
ylat = ((1:360)-181)/2; % grid cell latitude [ 90S-89.5N]
area = power(2*pi*earthRadius/(360*2),2)*cos(pi*ylat/180); % grid cell area
fnameout = rdir+'gridded_discharge_720x360x528.nc';  % output file name
%
nccreate(fnameout,'time','dimension',{'time',Inf},'Format','netcdf4');
ncwrite(fnameout,'time',tmn);
ncwriteatt(fnameout,'time','standard_name','time');
ncwriteatt(fnameout,'time','long_name','time');
ncwriteatt(fnameout,'time','units','months since 1979-12-15 00:00:00');
ncwriteatt(fnameout,'time','calendar','standard');
ncwriteatt(fnameout,'time','axis','T');
ncwriteatt(fnameout,'time','cell_methods','time: mean');
%
nccreate(fnameout,'lon','Dimensions',{'lon',720});
ncwrite(fnameout,'lon',xlon);
ncwriteatt(fnameout,'lon','standard_name','longitude');
ncwriteatt(fnameout,'lon','long_name','longitude');
ncwriteatt(fnameout,'lon','units','degrees_E');
ncwriteatt(fnameout,'lon','axis','X');
%
nccreate(fnameout,'lat','Dimensions',{'lat',360});
ncwrite(fnameout,'lat',ylat);
ncwriteatt(fnameout,'lat','standard_name','latitude');
ncwriteatt(fnameout,'lat','long_name','latitude');
ncwriteatt(fnameout,'lat','units','degrees_N');
ncwriteatt(fnameout,'lat','axis','Y');
%
nccreate(fnameout,'discharge','Dimensions',{'lon',720,'lat',360,'time',528},'FillValue',NaN);
ncwrite(fnameout,'discharge',discharge);
ncwriteatt(fnameout,'discharge','standard_name','discharge');
ncwriteatt(fnameout,'discharge','long_name','river_discharge');
ncwriteatt(fnameout,'discharge','units','m^3/s');
ncwriteatt(fnameout,'discharge','cell_methods','time: mean');
%
nccreate(fnameout,'area','Dimensions',{'lat',360});
ncwrite(fnameout,'area',area);
ncwriteatt(fnameout,'area','standard_name','grid_cell_area');
ncwriteatt(fnameout,'area','units','m^2');
%