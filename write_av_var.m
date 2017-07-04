% WRITE_BOX_AV
%
%  tims    Base [1990 1 1 +10], so tim = tims - greg2time([1990 1 1 10 0 0]);
%  bxid    If not just 0:(nbox-1)
%  fnm     output file name, incl '.nc'.  [default 'av_prop.nc']
%  vars    vector (ascending order):  1=temp  2=salt  3=w_flux  4=bx_nett
%                         5=eta
%  dat     cell array of matrices corresponding to 'vars'
%
% Jeff Dunn CMAR  19/1/06
%
% USAGE: write_box_av(tims,bxid,fnm,vars,dat)

function write_av_var(tims, bid, temp, salt, w, fnm)

nbx = length(bid);
ntm = length(tims);
nlv = size(temp,3);


% if min(tims) > 10000
%    t90  = greg2time([1990 1 1 10 0 0]);
%    tims = tims-t90;
% end

nc = netcdf(fnm,'noclobber');


% --- Create global attributes


% dimension variables
ncdim('time', 0, nc);
ncdim('level', nlv, nc);
ncdim('boxes', nbx, nc);
% time in seconds
nc{'time'}           = ncdouble('time');
nc{'time'}.long_name = 'time';
nc{'time'}.units     = 'days since 0001-01-01 00:00:00';
nc{'time'}.dt        = 1;
% polygons
nc{'boxes'}           = ncint('boxes');
nc{'boxes'}.long_name = 'Box IDs';
% layers
nc{'level'}           = ncint('level');
nc{'level'}.long_name = 'layer index; 1=near-surface';
nc{'level'}.positive  = 'down';


% Property variables
nc{'temperature'}             = ncfloat('time','boxes','level');
nc{'temperature'}.long_name   = 'temperature volume averaged';
nc{'temperature'}.units       = 'degree_C';
nc{'temperature'}.FillValue_  = -10e20;


nc{'salinity'}                = ncfloat('time','boxes','level');
nc{'salinity'}.long_name      = 'salinity volume averaged';
nc{'salinity'}.units          = '1e-3';
nc{'salinity'}.FillValue_     = -10e20;


nc{'verticalflux'}            = ncfloat('time','boxes','level');
nc{'verticalflux'}.long_name  = 'vertical flux averaged over floor of box';
nc{'verticalflux'}.positive   = 'upwards';
nc{'verticalflux'}.units      = 'm^3/s';
nc{'verticalflux'}.FillValue_ = -10e20;


%% General information
nc.Conventions = 'Standard';
nc.institution = 'CSIRO-IMAS';
nc.references = ' ';
str           = 'Properties averaged from hydrodynamic model outputs.';
nc.title      = str;
nc.source     = '';
nc.comment    = '';

% to mantain control and record
his        = ['Created by Bec Gorton and modified by Javier Porobic' date];
nc.history = his;

%% Data
for ii = 1:ntm
    nc{'temperature'}(ii,:,:)  = squeeze(temp(:, ii, :));
    nc{'salinity'}(ii,:,:)     = squeeze(salt(:, ii, :));
    nc{'verticalflux'}(ii,:,:) = squeeze(w(:, ii, :));
end

nc{'boxes'}(:) = bid;
nc{'level'}(:) = 1 : nlv;
nc{'time'}(:)  = tims;

%% close NETCDF file
nc             = close(nc);

%------------------------------------------------------------------------------
