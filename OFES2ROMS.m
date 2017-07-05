function OFES2ROMS(fnm, gfn)
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
nc     = netcdf.open(fnm);
%%ncg    = netcdf(grid);
lat    = netcdf.getVar(nc, netcdf.inqVarID(nc, 'lat'));
lon    = netcdf.getVar(nc, netcdf.inqVarID(nc,'lon')) - 360;
levels = 1 : length(netcdf.getVar(nc, netcdf.inqVarID(nc,'lev')));
s_rho  = tsmovavg(levels, 's', 2);
s_rho  = s_rho(~isnan(s_rho));
months = 1 : length(netcdf.getVar(nc, netcdf.inqVarID(nc,'time')));
v      = netcdf.getVar(nc, netcdf.inqVarID(nc,'v'));
u      = netcdf.getVar(nc, netcdf.inqVarID(nc,'u'));
w      = netcdf.getVar(nc, netcdf.inqVarID(nc,'w'));
sal    = netcdf.getVar(nc, netcdf.inqVarID(nc,'salinity'));
tem    = netcdf.getVar(nc, netcdf.inqVarID(nc,'temp'));
h      = netcdf.getVar(nc, netcdf.inqVarID(nc,'h'));
nc_t   = netcdf.getVar(nc, netcdf.inqVarID(nc,'time'));
nc_t   = nc_t  - 693961; % days since 1900 - 01 - 01
nc_t   = nc_t * 86400; % time in seconds
netcdf.close(fnm);
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
[lon, lat, lev, time] = ndgrid(lon, lat, levels, months) ;
[rho_lo, rho_la, rho_lev, rho_time] = ndgrid(rho_lon, rho_lat, s_rho, months) ;

%% from polar coordinates to Rho
v_rho    = interpn(lon, lat, lev, time, v, rho_lo, rho_la, rho_lev, rho_time);
u_rho    = interpn(lon, lat, lev, time, u, rho_lo, rho_la, rho_lev, rho_time);
w_rho    = interpn(lon, lat, lev, time, w, rho_lo, rho_la, rho_lev, rho_time);
temp_rho = interpn(lon, lat, lev, time, tem, rho_lo, rho_la, rho_lev, rho_time);
sal_rho  = interpn(lon, lat, lev, time, sal, rho_lo, rho_la, rho_lev, rho_time);
h        = interpn(squeeze(lon( : ,  : ,  1, 1)), squeeze(lat( : ,  : ,  1, 1)), h, squeeze(rho_lo( : ,  : ,  1, 1)), squeeze(rho_la( : ,  : ,  1, 1)));
mask     = h ; mask(mask < 0 | isnan(mask)) = 0; mask(mask > 0) = 1;

%% from cm/s to m/s
v_rho    = v_rho / 100;
u_rho    = u_rho / 100;
w_rho    = w_rho / 100;
%% creating a fake z - value
zeta     = squeeze(v_rho(:,  : , 1, :)) * 0;
zeta     = zeta  + 0.1;
%% removing NAs
v_rho(isnan(v_rho))       = 0;
u_rho(isnan(u_rho))       = 0;
w_rho(isnan(w_rho))       = 0;
temp_rho(isnan(temp_rho)) = -10e20;
sal_rho(isnan(sal_rho))   = -10e20;
h(h < 0 | isnan(h))       = 0;
%% Creating ncdf %%
nc_new = netcdf.create(gfn, '64BIT_OFFSET');
%% dimensions
xi_rho    =  netcdf.defDim(nc_new ,'xi_rho' , length(rho_lon));
eta_rho   =  netcdf.defDim(nc_new,'eta_rho', length(rho_lat));
time      =  netcdf.defDim(nc_new, 'time', size(time, 4));
levels    =  netcdf.defDim(nc_new, 'levels', size(lev, 3));
s_rho     =  netcdf.defDim(nc_new,'s_rho', size(s_rho, 2));
NC_GLOBAL = netcdf.getConstant('NC_GLOBAL'); % constant
%% the variables dimention is var[time, lev, eta, Xi]
%% VAR AND ATRIBUTES
var_scrum_time =  netcdf.defVar(nc_new, 'scrum_time', 'double', time);
var_mask       =  netcdf.defVar(nc_new, 'mask', 'float', [xi_rho eta_rho]);
var_h          =  netcdf.defVar(nc_new, 'h', 'float', [eta_rho  xi_rho]);
var_lat_rho    =  netcdf.defVar(nc_new, 'lat_rho', 'float', [xi_rho eta_rho]);
var_lon_rho    =  netcdf.defVar(nc_new, 'lon_rho', 'float', [xi_rho eta_rho]);
var_v          =  netcdf.defVar(nc_new, 'v', 'float', [xi_rho eta_rho s_rho time]);
var_u          =  netcdf.defVar(nc_new, 'u', 'float', [xi_rho eta_rho s_rho time]);
var_w          =  netcdf.defVar(nc_new, 'w', 'float', [xi_rho eta_rho s_rho time]);
var_salinity   =  netcdf.defVar(nc_new, 'salinity', 'float', [xi_rho eta_rho s_rho time]);
var_temp       =  netcdf.defVar(nc_new, 'temp', 'float', [xi_rho eta_rho s_rho time]);
var_zeta       =  netcdf.defVar(nc_new, 'zeta', 'float', [time xi_rho eta_rho]);
var_pm         =  netcdf.defVar(nc_new, 'pm', 'float', [xi_rho eta_rho]);
var_pn         =  netcdf.defVar(nc_new, 'pn', 'float', [xi_rho eta_rho]);
%netcdf.reDef(nc_new)
%% time_from
netcdf.putAtt(nc_new, var_scrum_time , 'long_name', char('scrum_time'));
netcdf.putAtt(nc_new, var_scrum_time , 'units', char('seconds since 1900-01-01 00:00:0.0'));
netcdf.putAtt(nc_new, var_scrum_time , 'grads_step', char('seconds'));
%% mask
netcdf.putAtt(nc_new, var_mask , 'long_name', char('mask on RHO-points'));
netcdf.putAtt(nc_new, var_mask , 'option_0', char('land'));
netcdf.putAtt(nc_new, var_mask , 'option_1', char('water'));

%% bathymetry
netcdf.putAtt(nc_new, var_h , 'long_name', char('bathymetry at RHO-points'));
netcdf.putAtt(nc_new, var_h , 'units', char('meter'));
netcdf.putAtt(nc_new, var_h , 'field', char('Bath,  scalar'));
%% Latitude
netcdf.putAtt(nc_new, var_lat_rho , 'long_name', char('lat_rho'));
netcdf.putAtt(nc_new, var_lat_rho , 'units', char('degree_north'));
netcdf.putAtt(nc_new, var_lat_rho , 'field', char('lat_rho,  scalar'));

%% longitude
netcdf.putAtt(nc_new, var_lon_rho , 'long_name', char('lon_rho'));
netcdf.putAtt(nc_new, var_lon_rho , 'units', char('degree_north'));
netcdf.putAtt(nc_new, var_lon_rho , 'field', char('lon_rho,  scalar'));

%% V - velocity component
netcdf.putAtt(nc_new, var_v , 'long_name', char('averaged v-momentum component'));
netcdf.putAtt(nc_new, var_v , 'units', char('meter second-1'));
netcdf.putAtt(nc_new, var_v , 'field', char('v-velocity,  scalar, series'));

%% U - velocity component
netcdf.putAtt(nc_new, var_u , 'long_name', char('averaged u-momentum component'));
netcdf.putAtt(nc_new, var_u , 'units', char('meter second-1'));
netcdf.putAtt(nc_new, var_u , 'field', char('v-velocity,  scalar, series'));

%% W - velocity component
netcdf.putAtt(nc_new, var_w , 'long_name', char('averaged v-momentum component'));
netcdf.putAtt(nc_new, var_w , 'units', char('meter second-1'));
netcdf.putAtt(nc_new, var_w , 'field', char('v-velocity,  scalar, series'));

%% Salinity
netcdf.putAtt(nc_new, var_salinity , 'long_name', char('salinity * 1000+35 [psu]'));
netcdf.putAtt(nc_new, var_salinity , 'units', char('practical salinity unit [psu]'));
netcdf.putAtt(nc_new, var_salinity , 'FillValue_', -10e20);

%% Temperature
netcdf.putAtt(nc_new, var_temp , 'long_name', char('potential temperature [c]'));
netcdf.putAtt(nc_new, var_temp , 'units', char('degrees celsius [c]'));
netcdf.putAtt(nc_new, var_temp , 'FillValue_', -10e20);

%% Zeta
netcdf.putAtt(nc_new, var_zeta , 'long_name', char('averaged free-surface [Fixed in 0.1 m]'));
netcdf.putAtt(nc_new, var_zeta , 'units', char('meter second-1'));
netcdf.putAtt(nc_new, var_zeta , 'field', char('free-surface,  scalar, series'));

%% pm
netcdf.putAtt(nc_new, var_pm , 'long_name', char('curvilinear coordinate metric in XI'));
netcdf.putAtt(nc_new, var_pm , 'units', char('meter second-1'));
netcdf.putAtt(nc_new, var_pm , 'field', char('pm, scalar'));

%% pn
netcdf.putAtt(nc_new, var_pn , 'long_name', char('curvilinear coordinate metric in ETA'));
netcdf.putAtt(nc_new, var_pn , 'units', char('meter second-1'));
netcdf.putAtt(nc_new, var_pn , 'field', char('pn, scalar'));

%% Global variables
str = ['ROMS File: Information extrapolated from the outputs of the OFES model ' ...
       '(http://www.jamstec.go.jp) using Javier Porobic approach (more detail : https://github.com/jporobicg)'];
his = ['Created Using Javier Porobic (more details : https://github.com/jporobicg) codes on' date];
netcdf.putAtt(nc_new, NC_GLOBAL, 'Title', char(' Roms files from OFES output'));
netcdf.putAtt(nc_new, NC_GLOBAL, 'Comments', str);
netcdf.putAtt(nc_new, NC_GLOBAL, 'history', str);
netcdf.putAtt(nc_new, NC_GLOBAL, 'hc', single(10));
netcdf.putAtt(nc_new, NC_GLOBAL, 'Cs_r', single(cs_r));
netcdf.putAtt(nc_new, NC_GLOBAL, 'Cs_w', single(cs_w));

%% Leaving define mode
netcdf.endDef(nc_new)

%% WRITE VARIABLES
netcdf.putVar(nc_new, var_scrum_time, nc_t);
netcdf.putVar(nc_new, var_mask, mask);
netcdf.putVar(nc_new, var_h, h);
netcdf.putVar(nc_new, var_lat_rho,  squeeze(rho_la(:, :, 1, 1)));
netcdf.putVar(nc_new, var_lon_rho,  squeeze(rho_lo(:, :, 1, 1)));
netcdf.putVar(nc_new, var_v, v_rho);
netcdf.putVar(nc_new, var_u, u_rho);
netcdf.putVar(nc_new, var_salinity, sal_rho);
netcdf.putVar(nc_new, var_temp, temp_rho);
netcdf.putVar(nc_new, var_zeta, zeta);
netcdf.putVar(nc_new, var_pm, pm);
netcdf.putVar(nc_new, var_pn, pn);
%% CLOSE FILES
netcdf.close(nc_new);