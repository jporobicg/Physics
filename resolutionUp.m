% Increase resoltuion of OFES hidrodynamic model for
%      JUAN FERNANDEZ RIDGE ECOSYSTEM MODEL
%
%%      Creadted by :  (Javier Porobic)  %%
%% OFES VERSION
% INPUT
% fnm  = Name of the OFES netcdf file
% grid = Name of the grid netcdf file for the ofes model
% sfn  = Name of the file to be saved
% USAGE: resolutionUp(fnm,grid);
function resolutionUp(fnm, grid)

%% reading variables
fil     =netcdf(fnm)
ncg     = netcdf(grid);
lat_ori = fil{'lat'}( : );
lon_ori = fil{'lon'}( : );
lev     = fil{'lev'}( : );
sal     = fil{'salinity'}(:, :);
tem     = fil{'temp'}(:, :);
u       = fil{'u'}(:, :);
w       = fil{'w'}(:, :);
v       = fil{'v'}(:, :);
time    = fil{'time'}(:, :);
h       = ncg{'ht'}(:);
%% Removing nans and missing values
v(v  == -999000000 | isnan(v)) = 0;
u(u  == -999000000 | isnan(u)) = 0;
w(w  == -999000000 | isnan(w)) = 0;
sal(sal  == -999000000 | isnan(sal)) = 0;
tem(tem  == -999000000 | isnan(tem)) = 0;
h(h  == -999000000 | isnan(h)) = 0;
%% increasing the resolution 5 times
lat_new = [min(lat_ori)  : 0.02  : max(lat_ori)];
lon_new = [min(lon_ori)  : 0.02  : max(lon_ori)];
%% removing land in selkirk and robinson
for t  = 1 : length(time)
    for l  = 1 : length(lev)
        %% selkirk
        tem(t, l, 53,  43 ) = mean(mean(tem(t, l, 52 : 2 : 54, 41 : 2 : 43)));
        sal(t, l, 53,  43 ) = mean(mean(sal(t, l, 52 : 2 : 54, 41 : 2 : 43)));
        v(t, l, 53,  43 )   = mean(mean(v(t, l, 52 : 2 : 54, 41 : 2 : 43)));
        u(t, l, 53,  43 )   = mean(mean(u(t, l, 52 : 2 : 54, 41 : 2 : 43)));
        w(t, l, 53,  43 )   = mean(mean(w(t, l, 52 : 2 : 54, 41 : 2 : 43)));
        %% Robinson
        tem(t, l, 54 ,  61 : 62) = mean(tem(t, l, [53 55] ,  61 : 62));
        sal(t, l, 54 ,  61 : 62) = mean(sal(t, l, [53 55] ,  61 : 62));
        v(t, l, 54 ,  61 : 62)   = mean(v(t, l, [53 55] ,  61 : 62));
        u(t, l, 54 ,  61 : 62)   = mean(u(t, l, [53 55] ,  61 : 62));
        w(t, l, 54 ,  61 : 62)   = mean(w(t, l, [53 55] ,  61 : 62));
    end
end
h(53,  43 )      = mean(mean(h(52 : 2 : 54, 41 : 2 : 43)));
h(54 ,  61 : 62) = mean(h([53 55] ,  61 : 62));
%% making grids
[time_gd, lev_gd, lat_gd, lon_gd]     = ndgrid(time, lev, lat_ori, lon_ori) ;
[lat_h, lon_h]                        = ndgrid(lat_ori, lon_ori) ;
[lat_h_new, lon_h_new]                = ndgrid(lat_new, lon_new) ;
[time_ngd, lev_ngd, lat_ngd, lon_ngd] = ndgrid(time, lev, lat_new, lon_new) ;
%% interpolation
v_ext    = interpn(time_gd, lev_gd, lat_gd, lon_gd, v, time_ngd, lev_ngd, lat_ngd, lon_ngd);
u_ext    = interpn(time_gd, lev_gd, lat_gd, lon_gd, u, time_ngd, lev_ngd, lat_ngd, lon_ngd);
w_ext    = interpn(time_gd, lev_gd, lat_gd, lon_gd, w, time_ngd, lev_ngd, lat_ngd, lon_ngd);
temp_ext = interpn(time_gd, lev_gd, lat_gd, lon_gd, tem, time_ngd, lev_ngd, lat_ngd, lon_ngd);
sal_ext  = interpn(time_gd, lev_gd, lat_gd, lon_gd, sal, time_ngd, lev_ngd, lat_ngd, lon_ngd);
h_ext    = interpn(lat_h, lon_h, h, lat_h_new, lon_h_new);
%% Creating NETCDF
nc_new = netcdf.create(gfn, '64BIT_OFFSET');
%% dimensions
lat_dim  = netcdf.defDim(nc_new , 'lat', size(lat_ngd, 3));
lev_dim  = netcdf.defDim(nc_new, 'lev', size(lev_gd, 2));
time_dim = netcdf.defDim(nc_new, 'time', size(time_gd, 1));
lon_dim  = netcdf.defDim(nc_new, 'lon', size(lon_ngd, 4));
%% Variables and attributes:
%% latitude
lat_var  = netcdf.defVar(nc_new, 'lat', 'double', lat_dim);
lon_var  = netcdf.defVar(nc_new, 'lon', 'double', lon_dim);
lev_var  = netcdf.defVar(nc_new, 'lev', 'double', lev_dim);
time_var = netcdf.defVar(nc_new, 'time', 'double', time_dim);
h_var    = netcdf.defVar(nc_new, 'h', 'float', [lon_dim lat_dim]);
temp_var = netcdf.defVar(nc_new, 'temp', 'float', [lon_dim lat_dim lev_dim time_dim]);
sal_var  = netcdf.defVar(nc_new, 'salinity', 'float', [lon_dim lat_dim lev_dim time_dim]);
u_var    = netcdf.defVar(nc_new, 'u', 'float', [lon_dim lat_dim lev_dim time_dim]);
v_var    = netcdf.defVar(nc_new, 'v', 'float', [lon_dim lat_dim lev_dim time_dim]);
w_var    = netcdf.defVar(nc_new, 'w', 'float', [lon_dim lat_dim lev_dim time_dim]);
netcdf.endDef(nc_new)
%% WRITE VARIABLES
netcdf.putVar(nc_new, time_var, time);
netcdf.putVar(nc_new, lon_var, lon_new);
netcdf.putVar(nc_new, lat_var, lat_new);
netcdf.putVar(nc_new, lev_var, lev);
netcdf.putVar(nc_new, h_var, permute(h_ext, [2 1]));
netcdf.putVar(nc_new, temp_var, permute(temp_ext, [4 3 2 1]));
netcdf.putVar(nc_new, sal_var, permute(sal_ext, [4 3 2 1]));
netcdf.putVar(nc_new, v_var, permute(v_ext, [4 3 2 1]));
netcdf.putVar(nc_new, u_var, permute(u_ext, [4 3 2 1]));
netcdf.putVar(nc_new, w_var, permute(w_ext, [4 3 2 1]));
netcdf.close(nc_new)
