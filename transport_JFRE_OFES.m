% TRANSPORT_JFRE  Calculate transports through defined faces, from
%      JUAN FERNANDEZ RIDGE ECOSYSTEM MODEL
%
%%      Modified by:  (Javier Porobic)  %%
%%      Based on Jeff Dunn code '2009'  %%
%% OFES VERSION
% INPUT
%   vert    Corners of the model
%   pt1    [nfc 2]  x,y of start point of each face
%   pt2    [nfc 2]  x,y of end point of each face
%   dlev   Number of the vertical levels in the model
%   lr     [nfc 2]  Id of box to left, ID of box to right, for each face (Not use
%           in this model)
%   dinc   [Optional]  face integration step length (km)  [default 0.1]  May
%               want to increase if all large boxes.
%   rimn   [Optional]  Min number of integration steps per face [default 3]
%
% OUTPUT
%  T       [nfc ntims ndep]  Sv, +ve to left from start
%  bxid    [nbx 1]  ID of boxes in bxnet {eg bxnet(N,:,:) refers to box bxid(N)}
%          NOTE: if lr incluudes box ID of -9, this is changed to max(bid)+1
%  bxnet   [nbx ntims ndep]  net transport for each box (hopefully near zero!)
%          +ve means gaining water.
%           Sv
%
% USAGE: [T,bxid,bxnet] = transport_curvi(pt1,pt2,lr,dinc,rimn);

function [T, tims] = transport_JFRE_OFES(vert, pt1, pt2, dlev, dinc, rimn, year, fnm)
%fnm, gnc);
%% Global Variables  %%
    nc   = netcdf(fnm);
    %gncf = netcdf(gnc);
    dint = diff(dlev);
    nlay = length(dint);
    tims = nc{'time'}(:);
    ntm  = length(tims);
    rimx = 400;
    if nargin<5 | isempty(rimn)
        rimn = 3;
    end
    % standard number of integration steps per face
    if nargin < 4 | isempty(dinc)
        dinc = 0.1;
    end
    file1 = (['1st_JFRE_first_Step.mat']);
    if ~exist(file1, 'file')  % This part need to be run only ones. its related
                              % with the general model configuration.
        disp(['Using hydro file ' fnm]);
        %% Mods
        gridDepth   = nc{'h'}(:, :);
        gridDepth(gridDepth == -999000000) = nan;
        rho_lo      = nc{'lon'}(:,:) - 360;
        rho_la      = nc{'lat'}(:,:);
        rho_lo      = repmat(rho_lo.', length(rho_la), 1);
        rho_la      = repmat(rho_la, 1, size(rho_lo, 2));
        sigmaValues = 1 : (length(nc{'lev'}(:,:)) - 1);       % this is the number of sigma-values layers

        %% this aproximation give me the approx depth at any sigma point  %%
        gridLayerDepths = zeros((length(sigmaValues) + 1), size(gridDepth, 1), size(gridDepth, 2));
        for i = 1 : size(gridDepth, 1)
            for j = 1 : size(gridDepth, 2)
                gridLayerDepths( :,  i, j) = [0 fliplr(sigma2zeta(gridDepth(i, j), 20, 5 , 0.6, length(sigmaValues)) ...
                                                       *  - 1)];
            end
        end
        layerValues = nan(nlay, length(sigmaValues), size(gridDepth, 1), size(gridDepth, 2));
        for layer = 1 : nlay
            for i = 1 : size(gridDepth, 1)
                for j = 1 : size(gridDepth, 2)
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
        disp(['using exiting data - ', file1])
        load(file1)
    end

    %% Polygons Faces %%
    ndps  = length(sigmaValues);
    nfc   = size(pt1, 1);
    T     = repmat(nan, [nfc ntm nlay]);
    file2=(['2nd_JFRE_second_Step.mat']);
    if ~exist(file2, 'file')  % This part need to be run only ones. its related
                              % with the general model configuration.
        disp('Creating new file JFRE_Second_Step.mat')
        % Prepare the integration grid for each face
        ngrd = 0;
        for ifc = 1 : nfc  % Number of faces
            x  = [pt1(ifc,1) pt2(ifc,1)];
            y  = [pt1(ifc,2) pt2(ifc,2)];
            yd = y(2) - y(1); %% Distance between points
            xd = x(2) - x(1);
            if abs(xd) > .00001
                lcor = cosd(mean(y));
                xd   = lcor * xd;
            else
                lcor = 1;
            end
            rdist     = 111.191 * sqrt(yd .^ 2 + xd .^ 2);  % distance in kilometres
            dirn(ifc) = cart2pol(xd, yd);                    % cartesian to polar (cilindrical) coordinates gives th angle
            rinc      = ceil(rdist / dinc);                   % To cope with Large domain:
            if rinc < rimn
                rinc = rimn;
            end
            if rinc > rimx
                rinc = rimx;
            end
            % disp(['face ' num2str(ifc) ' rinc = ' num2str(rinc)])
            % increase resolution in the faces
            yi   = yd / rinc;
            Y    = (y(1) + (yi / 2)) : yi : (y(2) - (yi / 2));
            xi   = xd / (lcor * rinc);
            X    = (x(1) + (xi / 2)) : xi : (x(2) - (xi / 2));
            ninc = length(Y);
            if ninc==0
                ninc = length(X);
                Y = repmat(mean(y), size(X));
            else
                if isempty(X)
                    X = repmat(mean(x), size(Y));
                elseif ninc > 2
                    % Adjustment for lat correction variation along face
                    lenX = X(end) - X(1);
                    xi   = xd / rinc;
                    xn   = xi / cosd(Y(1));
                    X(1) = x(1) + xn / 2;
                    for gg = 2 : ninc
                        X(gg) = X(gg-1) + xi / cosd(Y(gg));
                    end
                    % Correcting the last X
                    lenXnew = X(end) - X(1);
                    tmp     = X - X(1);
                    X       = X(1) + tmp * lenX / lenXnew;
                end
            end
            % savind result in a list
            ii       = ngrd + (1:ninc);
            ngrd     = ngrd + ninc;
            idx{ifc} = ii;
            flo(ii)  = X;
            fla(ii)  = Y;
            % In this data all levels are above the bottom.
            abovebot(:, ii) = ones(ndps, ninc);
            % Need to take into account all levels.
            mxdp(ifc) = nfc;
            % Volume interval is distance (changed to m 10^6 so transport in Sv)
            % by depth intervals.
            vint{ifc} = dint .* rdist / 1000;
        end
        save(file2, 'vint', 'abovebot', 'idx', 'flo', 'fla', 'ngrd', 'dirn')
    else
        disp(['using exiting data - ', file2])
        load(file2)
    end

    %% Creation of the Final File %%
    file3 = ([num2str(year), '_JFRE_third_Step.mat']);
    if ~exist(file3, 'file')
        disp('Creating new file JFRE_third_Step.mat')
        %% read steps inside the document  %%
        for id = 1 : ntm
            u     = squeeze(nc{'u'} (id, :, :, :));
            v     = squeeze(nc{'v'} (id, :, :, :));
            %% from cm/s to m/s (OFES output)
            u     = u ./ 100;
            v     = v ./ 100;
            U     = zeros(nlay, ngrd);
            V     = zeros(nlay, ngrd);
            udata = squeeze(u);
            vdata = squeeze(v);
            %% U an V devided by layer,  based in the Sigma Value
            for layer = 1 : nlay
                for k = 1 : length(sigmaValues)
                    uData_new(layer, k, :, :) = squeeze(udata(k, :, :)) .* squeeze(layerValues(layer, k, :, :));
                    vData_new(layer, k, :, :) = squeeze(vdata(k, :, :)) .* squeeze(layerValues(layer, k, :, :));
                end
            end
            final_uData = squeeze(nanmean(uData_new, 2));
            final_vData = squeeze(nanmean(vData_new, 2));
            for layer = 1 : nlay
                layer_u_data = squeeze(final_uData(layer, :, :));
                layer_v_data = squeeze(final_vData(layer, :, :));
                %% removind the Nans from U and V
                %U
                valid_mod_u  = ~isnan(layer_u_data);
                valid_u_lo   = rho_lo(valid_mod_u);
                valid_u_la   = rho_la(valid_mod_u);
                valid_u_data =  layer_u_data(valid_mod_u);
                %% creating the array (griddata)
                U(layer, :)  = griddata(valid_u_lo, valid_u_la, valid_u_data, flo, fla);
                %V
                valid_mod_v  = ~isnan(layer_v_data);
                valid_v_lo   = rho_lo(valid_mod_v);
                valid_v_la   = rho_la(valid_mod_v);
                valid_v_data =  layer_v_data(valid_mod_v);
                %% creating the array (griddata)
                V(layer, :) = griddata(valid_v_lo, valid_v_la, valid_v_data, flo, fla);
            end
            [dir2,spd] = cart2pol(U,V);
            %% U and V along each face
            for jj = 1 : nfc
                ii = idx{jj};
                tt = spd(:, ii) .* sin(dir2(:, ii) - dirn(jj));
                % For each layer
                for lvs = 1:nlay
                    for index = 1:size(ii)
                        value(index) =  mean(tt(lvs, index));
                    end
                    T(jj, id, lvs) = vint{jj}(lvs) .* mean(value);
                end
            end
        end
        save(file3, 'T', 'tims')
    else
        disp(['Using previous', file3])
        load(file3)
    end
    fprintf('\r');
    disp('Done       ')
end
