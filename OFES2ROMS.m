function OFES2ROMS(fnm, gfn, grid)
%% Function to extrapolate outpust from the OFES model to ROMS
% The main purpose of this code is to translate from the OFES grid  to the ROMS
% grid. Also, this code passes the variables from the grid of OFES to the grid of
% roms
%
% fnm  = Names of the OFES file
% gfn  = Name of the ROMS files that will be created (can't write the file if exist)
% grid =  nc file that contain bathimety information
%
%  USAGE OFES2ROMS(fnm, gfn, grid)
%
% Creator:  Javier Porobic
% Date   :  28/jun/2017
% github :  https://github.com/jporobicg

%% Open File and variables
nc     = netcdf(fnm);
ncg    = netcdf(grid);
lat    = nc{'lat'}(:);
lon    = nc{'lon'}(:) - 360;
levels = 1 : length(nc{'lev'}(:));
s_rho  = tsmovavg(levels, 's', 2);
s_rho  = s_rho(~isnan(s_rho));
months = 1 : length(nc{'time'}(:));
v      = nc{'v'}(:, :);
u      = nc{'u'}(:, :);
w      = nc{'u'}(:, :);
sal    = nc{'salinity'}(:,:);
tem    = nc{'temp'}(:,:);
h      = ncg{'ht'}(:);
nc_t   = nc{'time'}(:);
nc_t   = nc_t  - 693961; % days since 1900 - 01 - 01
nc_t   = nc_t * 86400; % time in seconds
%% removing nan and fill values
v(v  == -999000000 | isnan(v)) = 0;
u(u  == -999000000 | isnan(u)) = 0;
w(w  == -999000000 | isnan(w)) = 0;
sal(sal  == -999000000 | isnan(sal)) = 0;
tem(tem  == -999000000 | isnan(tem)) = 0;
%% creating new variable and defining dimentions
rho_lat = [min(lat) - 0.05; av2(lat); max(lat) + 0.05];
rho_lon = [min(lon) - 0.05; av2(lon); max(lon) + 0.05];
%% Adding curvilinear coordinate (1/distance(pt1,pt2)
S = referenceSphere('earth', 'm'); %% Oblate ellipsoid
for pos = 2 : length(rho_lat)
    pm(pos - 1) = 1 / (distance([rho_lat(pos - 1), rho_lon(1)], [rho_lat(pos), rho_lon(1) ...
                   ], S));
end
for pos = 2 : length(rho_lon)
    pn(pos - 1) = 1 / (distance([rho_lat(1), rho_lon(pos - 1)], [rho_lat(1), rho_lon(pos) ...
                   ], S));
end
pmt      = [pm(1), pm];
pnt      = [pn(1), pn];
[pn, pm] = meshgrid(pnt, pmt);
%% S - coordinate critical depths Those values are from the interpolation of the
%% roms grid (there is no information for the OFES model about this)
cs_w = [-1, -0.9278, -0.8274, -0.737, -0.6586, -0.5883, -0.5247, -0.4675, -0.4175, ...
        -0.373, -0.3328, -0.2965, -0.2647, -0.2365, -0.211, -0.1881, -0.1677, ...
        -0.1499, -0.1338,  -0.1192, -0.1063, -0.095, -0.0848, -0.0756, -0.0673, ...
        -0.0602, -0.0537, -0.0479, -0.0426, -0.038, -0.0339, -0.0302, -0.0269, ...
        -0.024, -0.0214, -0.019, -0.0168, -0.015, -0.0133, -0.0117, -0.0103, -0.0091, ...
        -0.008, -0.007, -0.006, -0.0052, -0.0044, -0.0037, -0.003, -0.0023, -0.0017, ...
        -0.0011, -0.0006, 0];
cs_r = tsmovavg(cs_w, 's', 2);
cs_r = cs_r(~isnan(cs_r));
%% new grids
[time, lev, lat, lon] = ndgrid(months, levels, lat, lon) ;
[rho_time, rho_lev, rho_lat, rho_lon] = ndgrid(months, s_rho, rho_lat, rho_lon) ;
%% from polar coordinates to Rho
v_rho    = interpn(time, lev, lat, lon, v, rho_time, rho_lev, rho_lat, rho_lon);
u_rho    = interpn(time, lev, lat, lon, u, rho_time, rho_lev, rho_lat, rho_lon);
w_rho    = interpn(time, lev, lat, lon, w, rho_time, rho_lev, rho_lat, rho_lon);
temp_rho = interpn(time, lev, lat, lon, tem, rho_time, rho_lev, rho_lat, rho_lon);
sal_rho  = interpn(time, lev, lat, lon, sal, rho_time, rho_lev, rho_lat, rho_lon);
h        = interpn(squeeze(lat(1, 1,  : ,  : )), squeeze(lon(1, 1,  : ,  : )), h, squeeze(rho_lat(1, 1,  : ,  : )), squeeze(rho_lon(1, 1,  : ,  : )));
mask     = h ; mask(mask < 0 | isnan(mask)) = 0; mask(mask > 0) = 1;

%% from cm/s to m/s
v_rho    = v_rho / 100;
u_rho    = u_rho / 100;
w_rho    = w_rho / 100;
%% creating a fake z - value
zeta     = squeeze(v_rho(:, 1, :, :)) * 0;
zeta     = zeta  + 0.1;
%% removing NAs
v_rho(isnan(v_rho))       = 0;
u_rho(isnan(u_rho))       = 0;
w_rho(isnan(w_rho))       = 0;
temp_rho(isnan(temp_rho)) = -10e20;
sal_rho(isnan(sal_rho))   = -10e20;
h(h < 0 | isnan(h))       = 0;
%% Creating ncdf %%
nc_new = netcdf(gfn, 'clobber');
%% dimensions
ncdim('xi_rho' , size(rho_lon, 4), nc_new);
ncdim('eta_rho', size(rho_lat, 3), nc_new);
ncdim('time', size(time, 1), nc_new);
ncdim('levels', size(lev, 2), nc_new);
ncdim('s_rho', size(s_rho, 2), nc_new);
%% the variables dimention is var[time, lev, eta, Xi]
%% VAR AND ATRIBUTES
%% time_from
nc_new{'scrum_time'}            = ncdouble('time');
nc_new{'scrum_time'}.long_name  = 'scrum_time';
nc_new{'scrum_time'}.units      = ncchar('seconds since 1900-01-01 00:00:0.0');
nc_new{'scrum_time'}.grads_step = ncchar('seconds');
%% mask
nc_new{'mask'}           = ncfloat('eta_rho', 'xi_rho');
nc_new{'mask'}.long_name = ncchar('mask on RHO-points');
nc_new{'mask'}.option_0  = ncchar('land');
nc_new{'mask'}.option_1  = ncchar('water');
%% bathymetry
nc_new{'h'}            = ncfloat('eta_rho', 'xi_rho');
nc_new{'h'}.long_name  = ncchar('bathymetry at RHO-points');
nc_new{'h'}.units      = ncchar('meter');
nc_new{'h'}.field      = ncchar('bath, scalar');
%% Latitude
nc_new{'lat_rho'}           = ncfloat('eta_rho', 'xi_rho');
nc_new{'lat_rho'}.long_name = ncchar('lat_rho');
nc_new{'lat_rho'}.units     = ncchar('degree_north');
nc_new{'lat_rho'}.field     = ncchar('lat_rho, scalar');
%% longitude
nc_new{'lon_rho'}           = ncfloat('eta_rho', 'xi_rho');
nc_new{'lon_rho'}.long_name = ncchar('lon_rho');
nc_new{'lon_rho'}.units     = ncchar('degree_north');
nc_new{'lon_rho'}.field     = ncchar('lon_rho, scalar');
%% V - velocity component
nc_new{'v'}           = ncfloat('time', 's_rho', 'eta_rho', 'xi_rho');
nc_new{'v'}.long_name = ncchar('averaged v-momentum component');
nc_new{'v'}.units     = ncchar('meter second-1');
nc_new{'v'}.field     = ncchar('v-velocity, scalar, series');
%% U - velocity component
nc_new{'u'}           = ncfloat('time', 's_rho', 'eta_rho', 'xi_rho');
nc_new{'u'}.long_name = ncchar('averaged u-momentum component');
nc_new{'u'}.units     = ncchar('meter second-1');
nc_new{'u'}.field     = ncchar('v-velocity, scalar, series');
%% W - velocity component
nc_new{'w'}           = ncfloat('time', 's_rho', 'eta_rho', 'xi_rho');
nc_new{'w'}.long_name = ncchar('averaged vertically v-momentum component');
nc_new{'w'}.units     = ncchar('meter second-1');
nc_new{'w'}.field     = ncchar('w-velocity, scalar, series');
%% Salinity
nc_new{'salinity'}            = ncfloat('time', 's_rho', 'eta_rho', 'xi_rho');
nc_new{'salinity'}.long_name  = ncchar('salinity * 1000+35 [psu] ');
nc_new{'salinity'}.units      = ncchar('practical salinity unit [psu]');
nc_new{'salinity'}.FillValue_ = -10e20;
%% Temperature
nc_new{'temp'}            = ncfloat('time', 's_rho', 'eta_rho', 'xi_rho');
nc_new{'temp'}.long_name  = ncchar('potential temperature [c]');
nc_new{'temp'}.units      = ncchar('degrees celsius [c]');
nc_new{'temp'}.FillValue_ = -10e20;
%% Zeta
nc_new{'zeta'}           = ncfloat('time', 'eta_rho', 'xi_rho');
nc_new{'zeta'}.long_name = ncchar('averaged free-surface [Fixed in 0.1 m]');
nc_new{'zeta'}.units     = ncchar('meter second-1');
nc_new{'zeta'}.field     = ncchar('free-surface, scalar, series');
%% pm
nc_new{'pm'}           = ncfloat('eta_rho', 'xi_rho');
nc_new{'pm'}.long_name = ncchar('curvilinear coordinate metric in XI');
nc_new{'pm'}.units     = ncchar('meter-1');
nc_new{'pm'}.field     = ncchar('pm, scalar');
%% pn
nc_new{'pn'}           = ncfloat('eta_rho', 'xi_rho');
nc_new{'pn'}.long_name = ncchar('curvilinear coordinate metric in ETA');
nc_new{'pn'}.units     = ncchar('meter-1');
nc_new{'pn'}.field     = ncchar('pn, scalar');
%% WRITE VARIABLES
nc_new{'scrum_time'}(:)  = nc_t;
nc_new{'mask'}(:, :)     = mask;
nc_new{'h'}(:, :)        = h;
nc_new{'lat_rho'}(:, :)  = squeeze(rho_lat(1, 1, :, :));
nc_new{'lon_rho'}(:, :)  = squeeze(rho_lon(1, 1, :, :));
nc_new{'v'}(:, :)        = v_rho;
nc_new{'u'}(:, :)        = u_rho;
nc_new{'w'}(:, :)        = w_rho;
nc_new{'salinity'}(:, :) = sal_rho;
nc_new{'temp'}(:, :)     = temp_rho;
nc_new{'zeta'}(:, :)     = zeta;
nc_new{'pm'}(:, :)       = pm;
nc_new{'pn'}(:, :)       = pn;
%% Global information
str = ['ROMS File: Information extrapolated from the outputs of the OFES model ' ...
       '(http://www.jamstec.go.jp) using Javier Porobic approach (more detail : https://github.com/jporobicg)'];
his = ['Created Using Javier Porobic (more detail : https://github.com/jporobicg) codes on' date];
nc_new.title   = str;
nc_new.history = his;
nc_new.hc      = ncfloat(10);
nc_new.Cs_r    = ncfloat(cs_r);
nc_new.Cs_w    = ncfloat(cs_w);
nc_new         = close(nc_new);