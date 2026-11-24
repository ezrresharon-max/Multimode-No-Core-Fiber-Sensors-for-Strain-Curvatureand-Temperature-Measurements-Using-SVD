%% Check cross-sensitivity

clearvars

ddir = fullfile("Data", "DATOS");

paramDir = {
    % Folder    File name pattern with parameter value                  % Measured parameter
    'CURVA977'  'sand2cm977nm a ([0-9.]+) cm_USB2F022031__(.+)\.txt'    'Curvature'
    'TEMP977'   'sand977 a ([0-9.]+) °C_USB2F022031__(.+)\.txt'         'Temperature'
    'STRAIN977' 'sand2cm977nm a ([0-9.]+) um_USB2F022031__(.+)\.txt'    'Strain'
    };

measFolders = {
    '1ra medicion'
    '2da medicion'
    '3ra medicion'
    '4ta medicion'
    '5ta medicion'
    };

measIdx = [1:5];
% Parameters array
param   = cell(numel(measIdx),size(paramDir,1));
% Spectral data
spdata  = cell(numel(measIdx),size(paramDir,1));
% Index to map data to parameter
dataIdx = cell(numel(measIdx),size(paramDir,1));
for k = 1:size(paramDir,1)
    for midx = 1:numel(measIdx)
        % Full path to folder
        pdir = fullfile(ddir, paramDir{k,1}, measFolders{measIdx(midx)});
        % Get files list
        pfiles = dir(pdir);
        % Remove directory entries
        pfiles([pfiles.isdir]) = [];
        % Remove dot files from the list
        pfiles(cellfun(@(c) c(1) == '.', {pfiles.name})) = [];
        % Extract strings with parameter values
        paramString = regexp({pfiles.name}, paramDir{k,2}, 'tokens');
        % Convert parameter strings to numbers
        param{midx, k} = cellfun(@(c) str2double(c{1}{1}), paramString);
        if k == 1 && midx == 1
            % Read wavelengths only from first file
            lambda = load(fullfile(pdir, pfiles(1).name));
            lambda(:,2) = [];
            % Select spectral range with sensor data
            selectLambda = lambda > 840 & lambda < 1000;
        end
        % Allocate space for spectral data
        spdata{midx, k} = zeros(sum(selectLambda), length(pfiles));
        % Sort files order by parameter value
        [param{midx, k}, sortIdx] = sort(param{midx, k});
        % Set data to parameter index
        dataIdx{midx, k} = k * ones(1,length(pfiles));
        for i = 1:length(pfiles)
            % Full path to data file
            fname = fullfile(pdir, pfiles(sortIdx(i)).name);
            % Read data
            dd = load(fname);
            % Make sure that wavelength array is the same as in first file
            assert(isequal(lambda, dd(:,1)), 'Wavelength data is incorrect in %s', fname);
            % Load spectral data
            spdata{midx, k}(:,i) = dd(selectLambda,2);
        end
    end
end

% Convert cell arrays to data matrices
spdata = cell2mat(reshape(spdata, 1, []));
dataIdx = cell2mat(reshape(dataIdx, 1, []));
param = cell2mat(reshape(param, 1, []));

% Get orthogonal spectral functions from full data set
[u,s,v] = svds(spdata, 4);
amp = v * s;
% Calculate normalized amplitudes
ampNorm = amp(:,2:end) ./ amp(:,1);

% Plot orthogonal spectral functions
figure(1),clf
p1 = plot(lambda(selectLambda), u);
p1(1).LineWidth = 3;
grid on
xlabel('Wavelength, nm')
legend(arrayfun(@(n) sprintf('U_{%d}', n), [1:size(u,2)], ...
    UniformOutput = false), Location='northwest')

dataIdxUniq = unique(dataIdx);

% Plot normalized amplitudes vs parameter
figure(555), clf
tl = tiledlayout(3,1);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';
% Placeholders for axes
axparam = gobjects(3,1);
% Degree of polynomial approximation
polydegr = 4;
% Number of points to visualize smooth polynomial approximation
polyN = 100;
% Arrays:
% Coefficients of polynomial approximation
polycoef = zeros(polydegr+1, size(ampNorm,2), numel(dataIdxUniq));
% Coefficients of derivatives approximation
dercoef = zeros(polydegr, size(ampNorm,2), numel(dataIdxUniq));
% Parameters
paramX = zeros(polyN, numel(dataIdxUniq));
% polynomial approximation of An/A1
polyappr = zeros(polyN, size(ampNorm,2), numel(dataIdxUniq));
% derivatives approximation of An/A1
derappr = zeros(polyN, size(ampNorm,2), numel(dataIdxUniq));
colors = colororder;
for k = 1:numel(dataIdxUniq)
    % Parameter ^ power to calculate polynomial coefficients
    basefun = param(dataIdx == dataIdxUniq(k))' .^ [0:polydegr];
    % polynomial coefficients
    polycoef(:,:,k) = basefun \ ampNorm(dataIdx == dataIdxUniq(k),:);
    % derivatives coefficients
    dercoef(:,:,k) = polycoef(2:end,:,k) .* [1:polydegr]';
    % calculation of smooth polynomial approximation and derivatives
    [parmin, parmax] = bounds(param(dataIdx == dataIdxUniq(k)));
    paramX(:,k) = linspace(parmin, parmax, polyN)';
    polyappr(:,:,k) = (paramX(:,k) .^ [0:polydegr]) * polycoef(:,:,k);
    derappr(:,:,k) = (paramX(:,k) .^ [0:polydegr-1]) * dercoef(:,:,k);
    % New axes
    axparam(k) = nexttile(tl);
    % Plot normalized amplitudes vs parameter (experimental)
    plot(axparam(k), param(dataIdx == dataIdxUniq(k)), ...
        ampNorm(dataIdx == dataIdxUniq(k),:), 'o')
    hold on
    set(gca, "ColorOrderIndex", 1)
    % Plot normalized amplitudes vs parameter (polynomial approximation)
    plot(axparam(k), paramX(:,k), polyappr(:,:,k))
    grid on
    xlabel(axparam(k), paramDir{dataIdxUniq(k),3})
    ylabel(axparam(k), 'A_n / A_1')
    legend(arrayfun(@(n) sprintf('A_{%d}/A_1', n), [2:size(u,2)], ...
        UniformOutput = false), Location='east')
end

% Normalized derivatives (represent directions of changes of An/A1)
derapprNorm = derappr ./ sqrt(sum(derappr .^ 2, 2));

% Test point at initial points
derapprNorm1 = squeeze(derapprNorm(1,:,:));
% Covariance matrix (cosines on angles between vectors of derivatives)
covm = derapprNorm1' * derapprNorm1;
disp(covm)

% 3D representation
plot3d = @(x,varargin) plot3(x(:,1),x(:,2),x(:,3),varargin{:});

figure(777), clf
for k = 1:numel(dataIdxUniq)
    plot3d(ampNorm(dataIdx == dataIdxUniq(k),:), '.')
    hold on
end
set(gca, "ColorOrderIndex", 1)
for k = 1:numel(dataIdxUniq)
    plot3d(polyappr(:,:,k))
end
set(gca, "ColorOrderIndex", 1)
for k = 1:numel(dataIdxUniq)
    paramMin = min(param(dataIdx == dataIdxUniq(k)));
    coord = ampNorm(dataIdx == dataIdxUniq(k) & param == paramMin,:);
    plot3d(coord, 'o')
    str = sprintf('  %s = %4.1f', paramDir{dataIdxUniq(k),3}, paramMin);
    text(coord(1,1), coord(1,2), coord(1,3), str)
end
% set(gca, "ColorOrderIndex", 1)
% for k = 1:numel(dataIdxUniq)
%     paramMax = max(param(dataIdx == dataIdxUniq(k)));
%     coord = ampNorm(dataIdx == dataIdxUniq(k) & param == paramMax,:);
%     plot3d(coord, 'o')
%     str = sprintf('  %s = %f', paramDir{dataIdxUniq(k),3}, paramMax);
%     text(coord(1,1), coord(1,2), coord(1,3), str)
% end
grid on
axis equal
