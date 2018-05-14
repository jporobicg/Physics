%Calculate the Average of variables from the CMIP5 model output
%%      Created by:  (Javier Porobic)     %%
%%                                        %%
%
% INPUT
%   vert   cell array, each cell has [n 2] vertices of a polygon
%   varn   'salt' or 'temp'
%
% OUTPUT
%  av      [nbox ntims ndep]
%  nav     [nbox]  number of model grid points per box
%
%
% USAGE: [av,nav] = box_av_curvi(vert,varn,tims);

function [av,nav] = box_avg_VAR_CC(vert, varn, dlev, year, fnm, gdirec)

nc   = netcdf(fnm);
tims = nc{'time'}(:);
dint = diff(dlev);
nlay = length(dint);
ntm  = length(tims);
nbox = length(vert);

file1=([gdirec, 'Var_CC_first_Step.mat']);
if ~exist(file1, 'file')  % This part need to be run only ones. its related
                          % with the general model configuration.
    disp(['Using hydro file ' fnm]);
    %gridDepth   = gncf{'ht'}(:, :);
    gridDepth   = nc{'h'}(:, :);
    %uniDepth    = squeeze(gridDepth(:,1,1));
    gridDepth(gridDepth == -999000000) = nan;
    rho_lo      = nc{'lon'}(:,:); %- 360;
    rho_la      = nc{'lat'}(:,:);
    rho_lo      = repmat(rho_lo.', length(rho_la), 1);
    rho_la      = repmat(rho_la, 1, size(rho_lo, 2));
    sigmaValues = 1 : (length(nc{'lev'}(:, :)) -1);

    %% this aproximation give me the approx depth at any sigma point  %%
    gridLayerDepths = gridDepth;
    layerValues     = nan(nlay, length(sigmaValues), size(gridDepth, 2), size(gridDepth, 3));
    for layer = 1 : nlay
        for i = 1 : size(gridDepth, 2)
            for j = 1 : size(gridDepth, 3)
                for k = 1 : length(sigmaValues)
                    minLayer = dlev(layer);
                    maxLayer = dlev(layer + 1);
                    minSigma = gridLayerDepths(k, i, j);
                    maxSigma = gridLayerDepths(k + 1, i, j);
                    if minLayer < maxSigma && maxLayer > minSigma
                        layerValues(layer, k, i, j) = 1;
                    end
                end
            end
        end
    end
    save(file1,  'gridLayerDepths', 'gridDepth', 'sigmaValues', 'rho_lo', 'rho_la', 'layerValues')
else
    disp(['Loading - ', file1])
    load(file1)
end

file2=([gdirec, varn, '_CC_Second_Step.mat']);
if ~exist(file2, 'file') % This part need to be run only ones. its related
                         % with the general model configuration.
    barea = zeros(nbox);
    ngrd  = 0;
    for ibx = 1 : nbox  % Prepare the integration grid for all faces
        x          = vert{ibx}(:, 1);
        y          = vert{ibx}(:, 2);
        ig         = find(inpolygon(rho_lo, rho_la, x, y)); % cell by polygon
        nav(ibx)   = length(ig);
        ii         = ngrd + (1 : nav(ibx));
        ngrd       = ngrd + nav(ibx);
        idx{ibx}   = ii;
        flo(ii)    = rho_lo(ig);
        fla(ii)    = rho_la(ig);
        [tix, tiy] = find(inpolygon(rho_lo,rho_la,x,y));
        % disp(['xi', num2str(length(tix))])
        ix{ibx}    = tix;
        iy{ibx}    = tiy;
        %% the area of the plygon is only important for production %%
        if strcmp(varn,'w')
            % this pas lat, lon to x, y.
            % in that way its easy and fast to calculate the area
            xx         = (x - mean(x)) .* latcor(y) * 111119;
            yy         = y * 111119;
            barea(ibx) = polyarea(xx, yy) ./ 1000000;
        else
            barea(ibx) = 1.0;
        end
    end
    save(file2, 'idx', 'flo', 'fla', 'nav', 'iy', 'ix', 'gridLayerDepths', ...
         'barea')
    disp(['Saving  - ', file2]);
else
    disp(['loading  - ', file2]);
    load(file2)
end

Var_avg = repmat(nan, [nbox ntm nlay]);
file3   = ([gdirec, varn, '_CC_third_Step.mat']);
if ~exist(file3, 'file') % This part need to be run only ones. its related
    for id = 1 : ntm
        varData      = squeeze(nc{varn} (id, :, :, :));
        time_varData = double(squeeze(varData));
        time_varData(time_varData <= -999000000) = NaN;
        for layer = 1 : nlay
            for k = 1 : length(sigmaValues)
                nVarD(layer, k, :, :) = squeeze(time_varData(k, :, :)) .* squeeze(layerValues(layer, k, :, :));
            end
        end
        % make NaN the zeros to calculate the mena
        nVarD(nVarD == 0) = NaN;
        FvarD             = squeeze(nanmean(nVarD, 2));
        FvarD             = double(FvarD);
        if strcmp(varn,'salinity')
            FvarD = FvarD .* 1000 + 35;
        end
        if strcmp(varn,'w')
            FvarD = FvarD ./ 100;
        end
        % Grid the Data and interpolate
        for layer = 1 : nlay
            interpolatedVarData(layer,:, :) = griddata(rho_lo, rho_la, squeeze(FvarD(layer, :, :)), flo, fla);
        end
        interpolatedVarData = squeeze(interpolatedVarData);
        for jj = 1 : nbox  %% By box
            tmp = zeros(nlay,nav(jj));
            for ij = 1 : nav(jj)  % For each grid cell in this box.
                for lvs = 1 : nlay  % by layer
                                    %Average by box by grid cell by layer
                    tmp(lvs, ij) = barea(jj) .* nanmean(interpolatedVarData(lvs, idx{jj}(ij)));
                end
            end
            for lvs = 1:nlay
                Var_avg(jj,id,lvs) = nanmean(tmp(lvs,:));
            end
        end
    end
    save(file3, 'Var_avg', 'tims')
    disp(['Saving  - ', file3]);
else
    disp(['loading  - ', file3]);
    load(file3)
end

%---------------------------------------------------------------------
