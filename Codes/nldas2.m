%% Extracting the data from any big NETCDF file

nc_file = ('/home/suhail/Downloads/cheswx_interp_prcp_19480101_20171231.nc');
%nc_file = uigetfile()


ncdisp(nc_file)

lat = ncread(nc_file,"lat");

lon = ncread(nc_file,"lon");

time = ncread(nc_file,"time");


% Finding the exact or nearby gridpoint for given lat and lon from the netcdf

latitude = 39.8729;

Longitude = -75.2437;


lat_index = find(lat>latitude,1,'first');

lon_index = find(lon>Longitude,1,'first');



% Extracting the data from the netcdf

prcp = ncread(nc_file,'prcp',[lon_index lat_index 1],[1 1 Inf]);

prcp_cheswx_v1 = squeeze(prcp(1,1,:));

%syntax goes like [Inf 1 1 ]= number of slicing of the dimensions

% Create the time dimension

%  t1 = datetime(1981,1,1,0,0,0);
%  t2 = datetime(2013,12,31,0,0,0);
%  t = t1:days(1):t2;
% 
% t = t(:);
% t= t(1:length(prcp));

%% Visualization

% plot(t,prcp)
% xlabel("Years",'FontSize',20)
% ylabel("Total Precipitation (kg/m^2)","FontSize",20)