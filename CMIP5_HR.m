% Increase resoltuion of OFES hidrodynamic model for
%      JUAN FERNANDEZ RIDGE ECOSYSTEM MODEL
%
%%      Creadted by :  (Javier Porobic)  %%
%% OFES VERSION
% INPUT
% fnm  = Name of the netcdf file with low resolution
% sfn  = Name of the file to be saved
%clear all
%fnm  = '/home/demiurgo/Documents/PhD/Oceanography/projections/Data_CMIP5/Proj_2011.nc'
%fnm2  = '/media/demiurgo/JPG_backup/Hydro_Atlantis/ROMS-OFES/JFRE_ROMS_1951.nc'
%sfn = '/media/demiurgo/JPG_backup/Hydro_Atlantis/Climate_Change_CMIP5/Proj_2011_EXTR.nc'

% USAGE: resolutionUp(fnm,grid);
function CMIP5_HR(fnm, sfn)
disp(['Increasing resoution OFES file [', fnm, ' ]'])
%% reading variables
fil  = netcdf.open(fnm, 'NOWRITE');
lat  = netcdf.getVar( fil, netcdf.inqVarID(fil, 'lat'));
lon  = netcdf.getVar( fil, netcdf.inqVarID(fil, 'lon'));
tem  = netcdf.getVar( fil, netcdf.inqVarID(fil, 'thetao'));
time = netcdf.getVar( fil, netcdf.inqVarID(fil, 'time'));
h    = netcdf.getVar( fil, netcdf.inqVarID(fil, 'lev'));
netcdf.close(fil);
%% REmoving NaNs and Low values this is important for the extraplation
tem(tem  >=  1e+20 | isnan(tem)) = nan;
%% From Kelvins to Celcius
tem = tem - 273.15;
%% Deffine dimentions
depth   = 1 : 50;
lon_vec = double(lon(:, 1));
lat_vec = double(lat(1, :));

%% Tempral interpolation,  to get the values per day
for i = 1 : (length(time) - 1)
    tstp   = time(i) + (0 : sum((diff(time(i : (i + 1))) - 1)));
    tim    = time(i : (i + 1));
    [lon_gd, lat_gd, depth_gd, time_gd]      = ndgrid(lon_vec, lat_vec, depth, tim);
    [lon_new, lat_new, depth_new, time_new]  = ndgrid(lon_vec, lat_vec, depth, tstp);
    tem_n    = tem( : ,  : ,  : , i : (i + 1));
    temp_ext = interpn(lon_gd, lat_gd, depth_gd, time_gd, tem_n, lon_new, lat_new, ...
                       depth_new, time_new);
    if i == 1
        tmp_fin = temp_ext;
        t_fin = tstp;
    else
        tmp_fin = cat(4, tmp_fin, temp_ext);
        t_fin = [t_fin, tstp];
    end
        clear tem_n temp_ext

end

%% increasing the resolution 5 times
lat_new = [min(lat_vec)  : 0.2  : max(lat_vec)];
lon_new = [min(lon_vec)  : 0.2  : max(lon_vec)];

%% making grids
[lon_gd, lat_gd, lev_gd, time_gd]     = ndgrid(lon_vec, lat_vec, depth, t_fin);
[lon_h_new, lat_h_new, h_ext]         = ndgrid(lon_new, lat_new, h) ;
[lon_ngd, lat_ngd, lev_ngd, time_ngd] = ndgrid(lon_new, lat_new, depth, t_fin) ;
%% interpolation
temp_ext = interpn(lon_gd, lat_gd, lev_gd, time_gd, tmp_fin, lon_ngd, lat_ngd, lev_ngd, time_ngd);
%h_ext    = interpn(lon_h, lat_h, h_ext, lon_h_new, lat_h_new);
%% Creating NETCDF
nc_new = netcdf.create(sfn, '64BIT_OFFSET');
%% dimensions
lat_dim  = netcdf.defDim(nc_new , 'lat', size(lat_ngd, 2));
lev_dim  = netcdf.defDim(nc_new, 'lev', size(lev_gd, 3));
time_dim = netcdf.defDim(nc_new, 'time', size(time_gd, 4));
lon_dim  = netcdf.defDim(nc_new, 'lon', size(lon_ngd, 1));
%% Variables and attributes:
%% latitude
lat_var  = netcdf.defVar(nc_new, 'lat', 'double', lat_dim);
lon_var  = netcdf.defVar(nc_new, 'lon', 'double', lon_dim);
lev_var  = netcdf.defVar(nc_new, 'lev', 'double', lev_dim);
time_var = netcdf.defVar(nc_new, 'time', 'double', time_dim);
h_var    = netcdf.defVar(nc_new, 'h', 'float', [lon_dim lat_dim lev_dim]);
temp_var = netcdf.defVar(nc_new, 'temp', 'float', [lon_dim lat_dim lev_dim time_dim]);
%sal_var  = netcdf.defVar(nc_new, 'salinity', 'float', [lon_dim lat_dim lev_dim time_dim]);
%u_var    = netcdf.defVar(nc_new, 'u', 'float', [lon_dim lat_dim lev_dim time_dim]);
%v_var    = netcdf.defVar(nc_new, 'v', 'float', [lon_dim lat_dim lev_dim time_dim]);
%w_var    = netcdf.defVar(nc_new, 'w', 'float', [lon_dim lat_dim lev_dim time_dim]);
netcdf.endDef(nc_new)
%% WRITE VARIABLES
netcdf.putVar(nc_new, time_var, t_fin);
netcdf.putVar(nc_new, lon_var, lon_new);
netcdf.putVar(nc_new, lat_var, lat_new);
netcdf.putVar(nc_new, lev_var, depth);
netcdf.putVar(nc_new, h_var, h_ext);
netcdf.putVar(nc_new, temp_var, temp_ext);
%netcdf.putVar(nc_new, sal_var, permute(sal_ext, [4 3 2 1]));
%netcdf.putVar(nc_new, v_var, permute(v_ext, [4 3 2 1]));
%netcdf.putVar(nc_new, u_var, permute(u_ext, [4 3 2 1]));
%netcdf.putVar(nc_new, w_var, permute(w_ext, [4 3 2 1]));
netcdf.close(nc_new)

disp(['Done!...Saving New OFES file [', sfn, ' ]'])