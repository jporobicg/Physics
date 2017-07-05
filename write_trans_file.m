%  tims    Base [1990 1 1 +10], so tim = tims - greg2time([1990 1 1 10 0 0]);
%
%  fcid    If faces not just 0:(nfc-1), eg if indexed by irealfaces, then
%          fcid = (irealfaces-1)
%  fnm     output file name, incl '.nc'.  [default 'trans.nc']
%
% USAGE: write_trans_file(pt1,pt2,lr,tims,T,fcid,fnm)

function write_trans_file(pt1,pt2,lr,tims,T,fcid,fnm)

nfc = size(pt1,1);

if size(T,1) ~= nfc
   error
end

%disp('Check "time" units are appropriate')

%tims
% length(tims);
% if min(tims)>10000
%     %disp('*** Converting from time0 of 1900 to 1990')
%   t90 = greg2time([1990 1 1 10 0 0]);
%   tims = tims-t90;
% end

ntm = length(tims);
if size(T,2) ~= ntm
   size(T,2)
    ntm
   error
end
nlv = size(T,3);

if nargin<6 | isempty(fcid)
   fcid = 0:(nfc-1);
end

if nargin<7 | isempty(fnm)
   fnm = 'trans.nc';
end

%nc = netcdf.create(fnm, 'noclobber');
%nc = close(nc);
nc = netcdf(fnm,'noclobber');

%nc = netcdf(fnm, 'write')
%if isempty(nc)
%  error('Cant open file - does it already exist?');
%end

% --- Create global attributes


% dimension variables

% Define dimensions:
ncdim('time', 0, nc);
ncdim('level', nlv, nc);
ncdim('faces', nfc, nc);


%nc('time') = 0;
%nc('level') = nlv;
%nc('faces') = nfc;

nc{'time'} = ncdouble('time');
nc{'time'}.long_name = 'time';
nc{'time'}.units = 'days since 0001-01-01 00:00:00';
nc{'time'}.dt = 1;

nc{'faces'} = ncint('faces');
nc{'faces'}.long_name = 'Face IDs';

nc{'level'} = ncint('level');
nc{'level'}.long_name = 'layer index; 1=near-surface';
nc{'level'}.positive = 'down';

% other reference variables
nc{'pt1_x'} = ncfloat('faces');
nc{'pt1_x'}.long_name = 'x coordinate of point 1 of face';
nc{'pt1_x'}.units = 'degree_east';

nc{'pt1_y'} = ncfloat('faces');
nc{'pt1_y'}.long_name = 'y coordinate of point 1 of face';
nc{'pt1_y'}.units = 'degree_north';

nc{'pt2_x'} = ncfloat('faces');
nc{'pt2_x'}.long_name = 'x coordinate of point 2 of face';
nc{'pt2_x'}.units = 'degree_east';

nc{'pt2_y'} = ncfloat('faces');
nc{'pt2_y'}.long_name = 'y coordinate of point 2 of face';
nc{'pt2_y'}.units = 'degree_north';

nc{'dest_boxid'} = ncint('faces');
nc{'dest_boxid'}.long_name = 'ID of destination box';
nc{'dest_boxid'}.units = 'id';

nc{'source_boxid'} = ncint('faces');
nc{'source_boxid'}.long_name = 'ID of source box';
nc{'source_boxid'}.units = 'id';


% Transports
nc{'transport'} = ncfloat('time','faces','level');
nc{'transport'}.long_name = 'flux across face';
%nc{'transport'}.units = 'Sverdrup (10^6 m^3/s)';
nc{'transport'}.units = 'm^3/s';
nc{'transport'}.comment = '+ve is to left, viewing from pt1 to pt2';
nc{'transport'}.FillValue_ = -10e20;

nc.Conventions = 'Standard';
nc.institution = 'CSIRO - IMAS';
nc.references = ' ';

str = ['Transports across faces of box model derived from hydrodynamic' ...
       ' model current field. Positive numbers are water gain in dest' ...
       ' box, loss in source box.'];
%disp(str);
%sad = input('Type any text to add to this "title" : ','s');
%if ~isempty(sad)
%   str = [str sad];
%end
nc.title = str;
nc.source = '';%input('Type "source" attribute: ','s');
nc.comment = ''; %input('Type "comment" attribute: ','s');

his = ['Created by Bec Gorton, modified by Javier Porobic on' date];
%disp(his)
%sad = input('Add to this "history" : ','s');
%if ~isempty(sad)
%   his = [his sad];
%end
nc.history = his;

for ii = 1:ntm
   nc{'transport'}(ii,:,:) = squeeze(T(:,ii,:));
end

nc{'faces'}(:) = fcid;
nc{'level'}(:) = 1:nlv;
nc{'pt1_x'}(:) = pt1(:,1);
nc{'pt1_y'}(:) = pt1(:,2);
nc{'pt2_x'}(:) = pt2(:,1);
nc{'pt2_y'}(:) = pt2(:,2);
nc{'dest_boxid'}(:) = lr(:,1);
nc{'source_boxid'}(:) = lr(:,2);
nc{'time'}(:) = tims;

nc = close(nc);
