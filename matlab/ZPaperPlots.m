%ZPAPERPLOTS Publication plots and all 24 diagnostics from SPOT 7.2.
% Save this file in SPOT/Visualization (or in the SPOT root folder).
% Run any single-acquisition example, then run ZPaperPlots.
% It reads Out.Configuration and Out.Exposures{k}; it never reruns SPOT.
%
% Ordinary examples, CSV pools and multiple exposures:
%   run(fullfile('Examples','RUN_Example_GroundSky.m'));
%   ZPaperPlots;
%
% Dataset example: it returns index, not Out. Set InputSource='dataset_index'
% below and choose DatasetRows=1 or 'all'. Acquisition MAT files must exist.
% Saved acquisition: set InputSource='mat' and InputMATFile to acquisition.mat.
% The explicit source choice prevents a dataset run from plotting an old Out.
%
% Defaults: three separate regime panels, all 24 diagnostics in the same
% 4x6 order as SPOT_PlotStages, and a 2x3 dataset figure. Export is 1200 dpi.
% Diagnostic figure size and typography are sized for the complete 4x6 grid.
% Dataset masks and boxes are read directly from the acquisition output.
%
% Outputs: PaperPlots reports source, exposure, saved files and figure handles.
% Only cfg and PaperPlots are assigned in the caller; C, E and Out are retained.
% Default save location: <acquisition folder>/PaperPlots.
% Existing paper files with the same name are replaced when saving is enabled.
% No emulator plot/save flags need to be enabled when Out is in memory.

%% USER SETTINGS: input and saving
cfg = struct();
cfg.ScriptRevision = 'SPOT-24-panels-fine-stars';
cfg.InputSource = 'workspace';           % 'workspace', 'mat', 'dataset_index'
cfg.InputMATFile = '';                    % used only for InputSource='mat'
cfg.DatasetRows = 1;                      % row(s) of index, or 'all'
cfg.ExposureIndex = 1;                    % exposure(s), or 'all'
cfg.SaveFigures = true;
cfg.SaveDir = '';                         % empty: acquisition folder/PaperPlots
cfg.SaveFormat = 'png';                   % 'png', 'tif', 'pdf', 'eps'
cfg.DPI = 1200;                         % raster export resolution
cfg.SaveFIG = false;                      % additionally save editable MATLAB FIG
cfg.ShowFigures = true;                   % false: export quietly, close afterwards
cfg.CloseExistingPaperFigures = true;     % closes only figures tagged SPOTPaperPlot

%% USER SETTINGS: which figures
cfg.RunSeparateRegimePanels = true;       % camera truth / detector / RA-DEC
cfg.RunOverview = false;                  % 1x2: camera truth + detector
cfg.RunFullRADECTruth = false;             % standalone full RA/DEC panel
cfg.RunDetailedStageDiagnostics = true;   % all 24 diagnostic panels, 4 rows x 6 columns
cfg.RunTargetInspection = false;          % 1x2: detector + target-region zoom
cfg.RunDatasetProducts = true;            % 2x3: detector, masks, boxes, summary

%% USER SETTINGS: publication figure geometry
cfg.SingleColWidthIn = 3.50;
cfg.DoubleColWidthIn = 7.16;
cfg.SingleFigureIn = [cfg.SingleColWidthIn 2.55];
cfg.TwoPanelFigureIn = [cfg.DoubleColWidthIn 3.10];
cfg.RegimePanelFigureIn = [2.35 2.10];
cfg.StageFigureIn = [cfg.DoubleColWidthIn 5.20]; % full 4x6 diagnostic grid
cfg.DatasetFigureIn = [cfg.DoubleColWidthIn 4.25];
cfg.ScreenScale = 120;                    % screen pixels per inch; export uses DPI

%% USER SETTINGS: typography (stage fonts sized for 24 panels)
cfg.FontName = 'Times New Roman';
cfg.AxisFontSize = 6.5;
cfg.LabelFontSize = 8;
cfg.TargetLabelFontSize = 6;
cfg.PanelLabelFontSize = 9.0;
cfg.StageAxisFontSize = 7.0;
cfg.StageLabelFontSize = 7.4;
cfg.StagePanelTitleFontSize = 7.6;
cfg.StageTargetLabelFontSize = 5.0;
cfg.DatasetAxisFontSize = 6.2;
cfg.DatasetLabelFontSize = 6.6;
cfg.DatasetPanelTitleFontSize = 7.0;

%% USER SETTINGS: titles, axes and display
cfg.TitlePrefix = '';
cfg.CaseName = '';                        % empty: actual geometric regime
cfg.ShowSingleAxisTitles = false;
cfg.ShowStagePanelTitles = true;
cfg.ShowStageMainTitle = false;
cfg.AddPanelLabels = true;
cfg.CompactSeparatePanels = false;
cfg.ApplyTightInset = true;
cfg.AxisLineWidth = 0.75;
cfg.UseLogForPositiveImages = true;        % optical/radiometric stage display only
cfg.ShowColorbars = false;
cfg.RADECFullScale = true;
cfg.RADECStride = 8;                      % subsample only the plotted FoV surface
cfg.ShowDetectorTruthOverlay = true;
cfg.DisplayPercentiles = [1 99.7];        % detector display, DN data unchanged
cfg.CameraStarSize_pt2 = 0.25;           % CAMERA grey-star marker area
cfg.CameraStarAlpha = 0.35;              % 0 = transparent, 1 = opaque
cfg.RADECStarSize_pt2 = 4;               % RA/DEC marker area; set 0.25 for finer dots
cfg.RADECStarAlpha = 0.5;
cfg.MaxCameraTruthStars = 8000;           % drawing limit only; Inf plots all
cfg.MaxRADECStarsFull = 8000;
cfg.StageColormap = gray(256);
cfg.TargetZoomPadding_pix = 80;
cfg.DatasetShowTextPanel = true;
cfg.DatasetMaxSummaryRows = 8;

%% RETRIEVE AND PLOT (normally leave unchanged)
switch lower(char(cfg.InputSource))
    case 'workspace'
        assert(exist('Out','var') == 1, 'SPOT:PaperInput', ...
            ['Run an example that produces Out first. After the dataset ', ...
             'example, select cfg.InputSource=''dataset_index''.']);
        PaperPlots = local_run(cfg, Out);
    case 'mat'
        PaperPlots = local_run(cfg, []);
    case 'dataset_index'
        assert(exist('index','var') == 1 && istable(index), ...
            'SPOT:PaperInput', 'Run the dataset example first to produce index.');
        PaperPlots = local_run(cfg, index);
    otherwise
        error('SPOT:PaperInput', 'Unknown cfg.InputSource: %s.', cfg.InputSource);
end

%% LOCAL FUNCTIONS
function report = local_run(cfg, input)
local_validate_settings(cfg);
source = lower(char(cfg.InputSource));
paths = {''};
if strcmp(source, 'mat')
    assert(strlength(string(cfg.InputMATFile)) > 0, 'SPOT:PaperInput', ...
        'Set cfg.InputMATFile to an exported acquisition.mat file.');
    paths = {char(cfg.InputMATFile)};
elseif strcmp(source, 'dataset_index')
    assert(ismember('directory', input.Properties.VariableNames), ...
        'SPOT:PaperInput', 'index must contain the directory column from SPOT.batch.');
    rows = local_indices(cfg.DatasetRows, height(input), 'DatasetRows');
    directories = string(input.directory(rows));
    assert(~any(ismissing(directories) | strlength(directories) == 0), ...
        'SPOT:PaperInput', 'Selected dataset rows have no saved results. Enable E.SaveResults and C.Output.SaveMAT.');
    paths = cellstr(fullfile(directories, 'acquisition.mat'));
end
% Validate inputs before closing previously generated paper figures.
if strcmp(source, 'workspace')
    local_validate_output(input);
else
    for n = 1:numel(paths)
        assert(isfile(paths{n}), 'SPOT:PaperInput', ...
            'Acquisition file not found: %s. Check saved paths and C.Output.SaveMAT.', paths{n});
    end
end
if cfg.CloseExistingPaperFigures
    delete(findall(groot, 'Type', 'figure', 'Tag', 'SPOTPaperPlot'));
end
report = struct('Source', {}, 'SampleID', {}, 'ExposureIndex', {}, ...
    'Regime', {}, 'SaveDirectory', {}, 'Files', {}, 'FigureHandles', {});
fprintf('RUN_SPOT_PaperPlots revision: %s\n', cfg.ScriptRevision);
for n = 1:numel(paths)
    if strcmp(source, 'workspace')
        acquisition = input;
    else
        loaded = load(paths{n}, 'Out');
        assert(isfield(loaded, 'Out'), 'SPOT:PaperInput', ...
            'The MAT file must contain the exported variable Out: %s.', paths{n});
        acquisition = loaded.Out;
        clear loaded
        local_validate_output(acquisition);
    end
    C = acquisition.Configuration;
    exposures = local_indices(cfg.ExposureIndex, numel(acquisition.Exposures), 'ExposureIndex');
    for k = exposures
        E = acquisition.Exposures{k};
        local_validate_exposure(C, E, cfg);
        plotCfg = cfg;
        plotCfg.SaveDir = local_save_directory(cfg, acquisition, paths{n});
        caseTag = local_case_tag(cfg, E);
        % Explicit shared output folders still distinguish dataset samples.
        if numel(paths) > 1 || strcmp(source, 'dataset_index')
            caseTag = sprintf('%s_%s_row%04d', caseTag, ...
                local_safe_tag(C.Dataset.SampleID), rows(n));
        end
        [files, handles] = local_plot_exposure(C, E, k, plotCfg, caseTag);
        r = struct('Source', source, 'SampleID', char(C.Dataset.SampleID), ...
            'ExposureIndex', k, 'Regime', E.Truth.GeometryRegime, ...
            'SaveDirectory', plotCfg.SaveDir, 'Files', {files}, 'FigureHandles', {handles});
        if ~isempty(paths{n}), r.Source = paths{n}; end
        report(end+1) = r; %#ok<AGROW>
        fprintf('Paper plots: %s, exposure %d, %d files saved.\n', caseTag, k, numel(files));
        if ~isempty(files), fprintf('  %s\n', plotCfg.SaveDir); end
    end
end
end

function local_validate_settings(cfg)
assert(ismember(lower(char(cfg.SaveFormat)), {'png','tif','tiff','pdf','eps'}), ...
    'SPOT:PaperSettings', 'SaveFormat must be png, tif, pdf or eps.');
validateattributes(cfg.DPI, {'numeric'}, {'scalar','integer','positive','finite'});
for name = {'CameraStarSize_pt2','RADECStarSize_pt2'}
    validateattributes(cfg.(name{1}), {'numeric'}, {'scalar','positive','finite'});
end
for name = {'CameraStarAlpha','RADECStarAlpha'}
    validateattributes(cfg.(name{1}), {'numeric'}, {'scalar','finite','>=',0,'<=',1});
end
validateattributes(cfg.RADECStride, {'numeric'}, {'scalar','integer','positive','finite'});
validateattributes(cfg.ScreenScale, {'numeric'}, {'scalar','positive','finite'});
validateattributes(cfg.TargetZoomPadding_pix, {'numeric'}, {'scalar','nonnegative','finite'});
validateattributes(cfg.DatasetMaxSummaryRows, {'numeric'}, {'scalar','integer','nonnegative','finite'});
for name = {'MaxCameraTruthStars','MaxRADECStarsFull'}
    v = cfg.(name{1});
    assert(isnumeric(v) && isscalar(v) && v >= 0 && (isinf(v) || v == fix(v)), ...
        'SPOT:PaperSettings', '%s must be a nonnegative integer or Inf.', name{1});
end
for name = {'SingleFigureIn','TwoPanelFigureIn','RegimePanelFigureIn','StageFigureIn','DatasetFigureIn'}
    validateattributes(cfg.(name{1}), {'numeric'}, {'vector','numel',2,'positive','finite'});
end
p = cfg.DisplayPercentiles;
assert(isnumeric(p) && numel(p) == 2 && all(isfinite(p)) && ...
    p(1) >= 0 && p(2) <= 100 && p(1) < p(2), 'SPOT:PaperSettings', ...
    'DisplayPercentiles must be increasing values between 0 and 100.');
end

function selected = local_indices(requested, count, name)
if (ischar(requested) || isstring(requested)) && strcmpi(requested, 'all')
    selected = 1:count;
else
    assert(isnumeric(requested) && isvector(requested) && ~isempty(requested) && ...
        all(isfinite(requested)) && all(requested == fix(requested)) && ...
        all(requested >= 1 & requested <= count), 'SPOT:PaperSelection', ...
        'cfg.%s must contain indices from 1 to %d, or ''all''.', name, count);
    selected = unique(reshape(requested, 1, []), 'stable');
end
assert(~isempty(selected), 'SPOT:PaperSelection', 'No items available for cfg.%s.', name);
end

function local_validate_output(Out)
assert(isstruct(Out) && isscalar(Out) && all(isfield(Out, {'Configuration','Exposures'})) && ...
    iscell(Out.Exposures) && ~isempty(Out.Exposures), 'SPOT:PaperInput', ...
    'Expected a new SPOT Out with Configuration and a nonempty Exposures cell array.');
end

function local_validate_exposure(C, E, cfg)
assert(all(isfield(E, {'Image_DN','Truth','Targets','StarTruth','CameraToGCRF', ...
    'TargetTimeUTC','MidUTC','StartUTC','EndUTC'})), 'SPOT:PaperInput', ...
    'Exposure is missing the image, geometry or target histories required for plotting.');
assert(isequal(size(E.Image_DN), [C.Camera.Ny C.Camera.Nx]), ...
    'SPOT:PaperInput', 'Detector image size does not match Out.Configuration.Camera.');
if cfg.RunDetailedStageDiagnostics
    assert(all(isfield(E, {'StarRate','SkyRate','EarthRate','TargetRate', ...
        'SceneRate','OpticalStages','DetectorStages'})) && ...
        all(isfield(E.Truth, {'Active','EarthMask','SkyMask','SolidAngle_sr', ...
        'GeometryRegime','Elevation_deg'})) && ...
        all(isfield(E.OpticalStages, {'Vignetting','Vignetted','PSF'})) && ...
        all(isfield(E.DetectorStages, {'ExpectedSignal_e','ExpectedDark_e', ...
        'ReadNoise_e','Analog_DN','SignalElectrons','BeforeFullWell_e', ...
        'FullWellSaturated','ADCSaturated'})), 'SPOT:PaperInput', ...
        'All 24 diagnostics require the complete acquisition Out, not PNG quicklooks.');
end
if cfg.RunDatasetProducts
    assert(isfield(E, 'Labels') && all(isfield(E.Labels, {'Binary','Instance','Semantic'})), ...
        'SPOT:PaperInput', 'Dataset panels require the recorded Labels masks.');
end
end

function folder = local_save_directory(cfg, Out, matFile)
if ~isempty(cfg.SaveDir)
    folder = char(cfg.SaveDir);
elseif ~isempty(matFile)
    folder = fullfile(fileparts(matFile), 'PaperPlots');
elseif isfield(Out, 'ExportDirectory') && ~isempty(Out.ExportDirectory)
    folder = fullfile(char(Out.ExportDirectory), 'PaperPlots');
else
    C = Out.Configuration;
    folder = fullfile(char(C.Output.Directory), char(C.Dataset.SampleID), 'PaperPlots');
end
end

function tag = local_case_tag(cfg, E)
if ~isempty(cfg.CaseName)
    tag = local_safe_tag(cfg.CaseName);
    return
end
switch E.Truth.GeometryRegime
    case 'ground_sky', tag = 'GroundSky';
    case 'space_sky', tag = 'SpaceSky';
    case 'space_limb', tag = 'EarthLimb';
    case 'space_earth', tag = 'EarthBackground';
    otherwise, tag = local_safe_tag(E.Truth.GeometryRegime);
end
end

function tag = local_safe_tag(value)
tag = regexprep(char(string(value)), '[^A-Za-z0-9_-]', '_');
if isempty(tag), tag = 'SPOT'; end
end

function [files, handles] = local_plot_exposure(C, E, k, cfg, tag)
files = {}; handles = gobjects(0,1);
if cfg.RunSeparateRegimePanels
    panel = cfg; panel.AddPanelLabels = false;
    panel.CompactNormalAxisLabels = cfg.CompactSeparatePanels;
    names = {'CameraTruth','Detector','RADEC'};
    for n = 1:3
        fig = local_new_figure(panel, names{n}, 'regimepanel');
        ax = axes('Parent', fig);
        switch n
            case 1, local_plot_camera_truth(ax, C, E, panel);
            case 2, local_plot_detector(ax, C, E, panel, true);
            case 3, local_plot_radec_truth(ax, C, E, panel);
        end
        [f, h] = local_finish(fig, cfg, sprintf('%s_%s_exp%02d', tag, names{n}, k));
        files = [files f]; handles = [handles; h]; %#ok<AGROW>
    end
end
if cfg.RunOverview
    fig = local_new_figure(cfg, 'Overview', 'twopanel');
    tl = tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');
    ax = nexttile(tl); local_plot_camera_truth(ax, C, E, cfg); local_panel_label(ax, '(a)', cfg);
    ax = nexttile(tl); local_plot_detector(ax, C, E, cfg, true); local_panel_label(ax, '(b)', cfg);
    [f, h] = local_finish(fig, cfg, sprintf('%s_overview_exp%02d', tag, k));
    files = [files f]; handles = [handles; h];
end
if cfg.RunFullRADECTruth
    fig = local_new_figure(cfg, 'Full RA/DEC truth', 'single');
    ax = axes('Parent', fig); local_plot_radec_truth(ax, C, E, cfg);
    [f, h] = local_finish(fig, cfg, sprintf('%s_RADEC_truth_exp%02d', tag, k));
    files = [files f]; handles = [handles; h];
end
if cfg.RunDetailedStageDiagnostics
    fig = local_new_figure(cfg, 'Stage diagnostics', 'stage');
    local_plot_stages(fig, C, E, cfg);
    [f, h] = local_finish(fig, cfg, sprintf('%s_stage_diagnostics_exp%02d', tag, k));
    files = [files f]; handles = [handles; h];
end
if cfg.RunTargetInspection
    fig = local_new_figure(cfg, 'Target inspection', 'twopanel');
    tl = tiledlayout(fig, 1, 2, 'TileSpacing','compact', 'Padding','compact');
    ax = nexttile(tl); local_plot_detector(ax, C, E, cfg, true); local_panel_label(ax, '(a)', cfg);
    ax = nexttile(tl); local_plot_target_zoom(ax, C, E, cfg); local_panel_label(ax, '(b)', cfg);
    [f, h] = local_finish(fig, cfg, sprintf('%s_target_inspection_exp%02d', tag, k));
    files = [files f]; handles = [handles; h];
end
if cfg.RunDatasetProducts
    fig = local_new_figure(cfg, 'Dataset products', 'dataset');
    local_plot_dataset(fig, C, E, cfg);
    [f, h] = local_finish(fig, cfg, sprintf('%s_dataset_products_exp%02d', tag, k));
    files = [files f]; handles = [handles; h];
end
end

function fig = local_new_figure(cfg, name, style)
switch style
    case 'regimepanel', sz = cfg.RegimePanelFigureIn;
    case 'twopanel', sz = cfg.TwoPanelFigureIn;
    case 'stage', sz = cfg.StageFigureIn;
    case 'dataset', sz = cfg.DatasetFigureIn;
    otherwise, sz = cfg.SingleFigureIn;
end
visible = 'off'; if cfg.ShowFigures, visible = 'on'; end
fig = figure('Name',['SPOT paper plot | ' name], 'Tag','SPOTPaperPlot', ...
    'Color','w', 'Visible',visible, 'Units','pixels', ...
    'Position',[100 100 sz*cfg.ScreenScale], 'PaperPositionMode','manual', ...
    'PaperUnits','inches', 'PaperSize',sz, 'PaperPosition',[0 0 sz], ...
    'InvertHardcopy','off');
end

function [files, handles] = local_finish(fig, cfg, name)
% print honours physical PaperPosition and DPI independently of ScreenScale.
files = {}; handles = gobjects(0,1);
if cfg.SaveFigures
    if ~isfolder(cfg.SaveDir)
        [ok, message] = mkdir(cfg.SaveDir);
        assert(ok, 'SPOT:PaperExport', 'Cannot create output folder: %s', message);
    end
    fmt = lower(char(cfg.SaveFormat));
    switch fmt
        case 'png', device = '-dpng';
        case {'tif','tiff'}, device = '-dtiff';
        case 'pdf', device = '-dpdf';
        case 'eps', device = '-depsc';
    end
    path = fullfile(cfg.SaveDir, [name '.' fmt]);
    drawnow;
    print(fig, path, device, sprintf('-r%d', cfg.DPI));
    files{end+1} = path;
    if cfg.SaveFIG
        path = fullfile(cfg.SaveDir, [name '.fig']);
        savefig(fig, path); files{end+1} = path;
    end
end
if cfg.ShowFigures, handles = fig; else, close(fig); end
end

function local_plot_camera_truth(ax, C, E, cfg)
D = E.Truth;
support = .07*double(D.Active) + .13*double(D.EarthMask);
imagesc(ax, D.Xi_deg(1,:), D.Eta_deg(:,1), support);
axis(ax, 'image'); set(ax, 'YDir','normal', 'Color',[.015 .015 .02]);
colormap(ax, gray(256)); clim(ax, [0 1]); hold(ax, 'on');
S = E.StarTruth;
q = local_star_indices(numel(S.Xi_deg), cfg.MaxCameraTruthStars);
scatter(ax, S.Xi_deg(q), S.Eta_deg(q), cfg.CameraStarSize_pt2, [.65 .65 .65], ...
    'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',cfg.CameraStarAlpha);
half = C.Camera.FOV_deg/2;
if strcmp(C.Camera.Shape, 'circular')
    th = linspace(0, 2*pi, 361);
    plot(ax, atand(tand(half(1))*cos(th)), atand(tand(half(1))*sin(th)), ...
        '-', 'Color',[0 .85 0], 'LineWidth',1.4);
else
    rectangle(ax, 'Position',[-half 2*half], 'EdgeColor',[0 .85 0], 'LineWidth',1.4);
end
plot(ax, 0, 0, 'wx', 'LineWidth',2, 'MarkerSize',8);
local_plot_targets(ax, C, E, cfg, 'camera', 0, true, false);
xlim(ax, [-half(1) half(1)]*1.05); ylim(ax, [-half(2) half(2)]*1.05);
xlabel(ax, 'Cross-boresight camera angle, \xi [deg]');
ylabel(ax, 'Vertical camera angle, \eta [deg]');
local_single_title(ax, cfg, E, 'Camera-angle truth');
local_style(ax, cfg, 'normal'); local_compact(ax, cfg);
hold(ax, 'off');
end

function local_plot_detector(ax, C, E, cfg, overlay)
imagesc(ax, E.Image_DN); axis(ax, 'image'); set(ax, 'YDir','reverse');
colormap(ax, gray(256)); clim(ax, local_limits(E.Image_DN, cfg.DisplayPercentiles));
if cfg.ShowColorbars
    cb = colorbar(ax); cb.Label.String = 'Detector signal [DN]';
end
if overlay && cfg.ShowDetectorTruthOverlay
    local_plot_targets(ax, C, E, cfg, 'detector', 0, true, false);
end
xlabel(ax, 'x pixel'); ylabel(ax, 'y pixel');
local_single_title(ax, cfg, E, 'Detector image plane');
local_style(ax, cfg, 'normal'); local_compact(ax, cfg);
end

function local_plot_radec_truth(ax, C, E, cfg)
% The footprint and stars are the actual midpoint truth retained by SPOT.
% FullScale expands the coordinate axes; it does not fabricate all-sky stars.
[~, mid] = min(abs(seconds(E.TargetTimeUTC - E.MidUTC)));
b = E.CameraToGCRF(:,3,mid);
ra0 = mod(atan2d(b(2), b(1)), 360); dec0 = asind(max(-1, min(1, b(3))));
D = E.Truth;
ix = unique([1:cfg.RADECStride:C.Camera.Nx C.Camera.Nx]);
iy = unique([1:cfg.RADECStride:C.Camera.Ny C.Camera.Ny]);
ra = local_wrap(D.RA_deg(iy,ix) - ra0); dec = D.DEC_deg(iy,ix);
shade = .08 + .13*double(D.EarthMask(iy,ix));
shade(~D.Active(iy,ix)) = NaN;
seam = [false(size(ra,1),1) abs(diff(ra,1,2)) > 180] | ...
    [false(1,size(ra,2)); abs(diff(ra,1,1)) > 180];
ra(seam) = NaN;
surface(ax, ra, dec, zeros(size(ra)), shade, 'EdgeColor','none', 'FaceColor','flat');
view(ax, 2); colormap(ax, gray(256)); clim(ax, [0 1]);
set(ax, 'YDir','normal', 'Color',[.015 .015 .025]); hold(ax, 'on');
S = E.StarTruth; q = local_star_indices(numel(S.RA_deg), cfg.MaxRADECStarsFull);
scatter(ax, local_wrap(S.RA_deg(q)-ra0), S.DEC_deg(q), cfg.RADECStarSize_pt2, [.7 .7 .7], ...
    'filled', 'MarkerEdgeColor','none', 'MarkerFaceAlpha',cfg.RADECStarAlpha);
local_plot_targets(ax, C, E, cfg, 'radec', ra0, true, false);
% Marker size control in RA/DEC 
plot(ax, 0, dec0, 'wx', 'LineWidth',2, 'MarkerSize',2);
if cfg.RADECFullScale
    xlim(ax, [-180 180]); ylim(ax, [-90 90]);
else
    xlim(ax, [-.7 .7]*C.Camera.FOV_deg(1));
    ylim(ax, [max(-90, dec0-.7*C.Camera.FOV_deg(2)), min(90, dec0+.7*C.Camera.FOV_deg(2))]);
end
xlabel(ax, sprintf('Right ascension, \\alpha - %.2f^{\\circ} [deg]', ra0));
ylabel(ax, 'Declination, \delta [deg]');
local_single_title(ax, cfg, E, 'RA/DEC truth');
local_style(ax, cfg, 'normal');
set(ax, 'GridColor',[.55 .55 .55], 'GridAlpha',.22); local_compact(ax, cfg);
hold(ax, 'off');
end

function local_plot_stages(fig, C, E, cfg)
% The 24 data expressions and their order match Visualization/SPOT_PlotStages.
% Geometry masks are at the exposure midpoint. Scene component rates are
% exposure averages; optical arrays and detector stages are integrated.
% SkyMask describes celestial support, not visibility of targets over Earth.
D = E.Truth; O = E.OpticalStages; Q = E.DetectorStages;
blocked = D.EarthMask;
if strcmp(D.GeometryRegime, 'ground_sky')
    blocked = D.Active & D.Elevation_deg < C.Visibility.MinElevation_deg;
end
sky = E.SkyRate ./ max(D.SolidAngle_sr, realmin);
earth = E.EarthRate ./ max(D.SolidAngle_sr, realmin);
bias = 0; if C.Detector.BiasEnabled, bias = C.Detector.Bias_DN; end
% Panel 19 removes the configured bias and expected dark signal. It does
% not reverse saturation, quantization or a particular dark-noise draw.
data = {D.Active, blocked, D.SkyMask, D.EarthMask, E.StarRate, sky, ...
    earth, E.TargetRate, E.StarRate+E.TargetRate, sky+earth, E.SceneRate, O.Vignetting, ...
    O.Vignetted, O.PSF, Q.ExpectedSignal_e, Q.ExpectedDark_e, Q.ReadNoise_e, Q.Analog_DN, ...
    E.Image_DN-bias-Q.ExpectedDark_e/C.Detector.ConversionGain_e_per_DN, E.Image_DN, ...
    Q.SignalElectrons, Q.BeforeFullWell_e, Q.FullWellSaturated|Q.ADCSaturated, zeros(size(E.Image_DN))};
titles = {'In-FoV mask', 'Earth / horizon blocked mask', 'Visible sky mask', ...
    'Earth-in-FoV mask', 'Stars flux only', 'Sky background radiance', ...
    'Earth radiance only', 'Targets flux only', 'Point-source flux input', ...
    'Diffuse radiance input', 'Radiometric input to optics', 'Vignetting map', ...
    'After vignetting', 'After optics / PSF', 'Electron mean image', ...
    'Dark current mean', 'Read noise realization', 'Pre-ADC signal', ...
    'Bias/dark-subtracted image', 'Final detector image', 'After shot noise', ...
    'After dark current', 'Saturation mask', 'Detector truth map'};
units = {'Active FoV [0/1]', 'Earth/horizon blockage [0/1]', 'Visible sky [0/1]', ...
    'Earth in FoV [0/1]', 'Star rate [model units/s]', 'Sky radiance [model units/s/sr]', ...
    'Earth radiance [model units/s/sr]', 'Target rate [model units/s]', ...
    'Point-source rate [model units/s]', 'Diffuse radiance [model units/s/sr]', ...
    'Scene rate [model units/s]', 'Relative transmission', ...
    'Integrated optical signal [model units]', 'Integrated optical signal [model units]', ...
    'Expected signal [electrons]', 'Expected dark signal [electrons]', ...
    'Read noise [electrons]', 'Pre-ADC signal [DN]', 'Bias/dark-subtracted signal [DN]', ...
    'Detector signal [DN]', 'Signal after shot noise [electrons]', ...
    'Signal plus dark electrons', 'Saturated pixels [0/1]', 'Target truth'};
tl = tiledlayout(fig, 4, 6, 'TileSpacing','compact', 'Padding','compact');
for j = 1:24
    ax = nexttile(tl, j);
    isLog = ismember(j, [5:11 13 14]) && cfg.UseLogForPositiveImages;
    local_plot_stage(ax, data{j}, titles{j}, units{j}, isLog, cfg);
    % Fixed binary limits distinguish an all-zero mask from an all-one mask.
    if ismember(j, [1:4 23 24]), clim(ax, [0 1]); end
    if j == 24
        truthCfg = cfg; truthCfg.TargetLabelFontSize = cfg.StageTargetLabelFontSize;
        count = local_plot_targets(ax, C, E, truthCfg, 'detector', 0, true, false);
        if count == 0
            text(ax, .5, .5, {'no illuminated'; 'visible target'}, 'Units','normalized', ...
                'HorizontalAlignment','center', 'Color','w', 'FontSize',cfg.StageTargetLabelFontSize);
        end
    end
end
if cfg.ShowStageMainTitle
    title(tl, [cfg.TitlePrefix 'Stage-wise signal-chain traceability'], ...
        'Interpreter','none', 'FontSize',11, 'FontWeight','bold');
end
end

function local_plot_stage(ax, values, name, units, isLog, cfg)
a = double(values); a(~isfinite(a)) = NaN;
if isLog
    pos = a(isfinite(a) & a > 0);
    if ~isempty(pos)
        % Same positive-image display transform as SPOT_PlotStages.
        scale = max(max(pos)*1e-5, realmin);
        a = log10(1 + max(0,a)/scale);
        units = ['log10(1 + signal/scale); ' units];
    end
end
imagesc(ax, a); axis(ax, 'image'); set(ax, 'YDir','reverse');
colormap(ax, cfg.StageColormap); clim(ax, local_limits(a, [0 100]));
if cfg.ShowStagePanelTitles
    title(ax, local_stage_title(name), 'Interpreter','none');
end
if cfg.ShowColorbars, cb = colorbar(ax); cb.Label.String = units; end
xlabel(ax, 'x pixel'); ylabel(ax, 'y pixel'); local_style(ax, cfg, 'stage');
end

function lines = local_stage_title(name)
% Wrap long titles to fit six columns at the configured publication width.
words = strsplit(name); lines = {words{1}};
for j = 2:numel(words)
    if numel(lines{end}) + 1 + numel(words{j}) <= 23
        lines{end} = [lines{end} ' ' words{j}];
    else
        lines{end+1} = words{j}; %#ok<AGROW>
    end
end
end

function local_plot_target_zoom(ax, C, E, cfg)
xAll = []; yAll = [];
for j = 1:numel(E.Targets)
    [x, y, valid] = local_target_coordinates(E.Targets(j), 'detector', 0);
    xAll = [xAll x(valid)]; yAll = [yAll y(valid)]; %#ok<AGROW>
end
local_plot_detector(ax, C, E, cfg, true);
if isempty(xAll)
    local_single_title(ax, cfg, E, 'No visible target available for zoom');
    return
end
pad = max(1, cfg.TargetZoomPadding_pix);
xlim(ax, [max(.5, min(xAll)-pad), min(C.Camera.Nx+.5, max(xAll)+pad)]);
ylim(ax, [max(.5, min(yAll)-pad), min(C.Camera.Ny+.5, max(yAll)+pad)]);
local_single_title(ax, cfg, E, 'Zoomed target view');
end

function local_plot_dataset(fig, C, E, cfg)
% Exact masks: no dilation, threshold changes or independent label generator.
tl = tiledlayout(fig, 2, 3, 'TileSpacing','compact', 'Padding','compact');
ax = nexttile(tl, 1); local_plot_detector(ax, C, E, cfg, true);
title(ax, 'Detector image', 'Interpreter','none'); local_style(ax, cfg, 'dataset');
ax = nexttile(tl, 2);
local_plot_mask(ax, E.Labels.Binary, gray(256), [0 1], 'Binary target mask', cfg);
% Compact colour indices preserve arbitrary target IDs in the stored mask.
ids = unique(E.Labels.Instance(:)); ids(ids == 0) = [];
im = zeros(size(E.Labels.Instance));
for j = 1:numel(ids), im(E.Labels.Instance == ids(j)) = j; end
palette = [0 0 0; parula(max(1,numel(ids)))];
ax = nexttile(tl, 3);
local_plot_mask(ax, im, palette, [0 max(1,numel(ids))], 'Instance mask', cfg);
ax = nexttile(tl, 4);
local_plot_mask(ax, E.Labels.Semantic, gray(256), [0 2], 'Semantic mask: object / trail', cfg);
ax = nexttile(tl, 5); local_plot_detector(ax, C, E, cfg, false);
local_plot_targets(ax, C, E, cfg, 'detector', 0, false, true);
local_draw_boxes(ax, C, E, cfg);
title(ax, 'Detector, truth and bounding boxes', 'Interpreter','none');
local_style(ax, cfg, 'dataset');
ax = nexttile(tl, 6); axis(ax, 'off');
if cfg.DatasetShowTextPanel
    text(ax, .02, .98, local_summary(C, E, cfg), 'Units','normalized', ...
        'HorizontalAlignment','left', 'VerticalAlignment','top', ...
        'FontName','FixedWidth', 'FontSize',cfg.DatasetLabelFontSize, 'Interpreter','none');
end
title(ax, 'Label and metadata summary', 'Interpreter','none');
local_style(ax, cfg, 'dataset');
end

function local_plot_mask(ax, data, palette, limits, name, cfg)
imagesc(ax, data); axis(ax, 'image'); set(ax, 'YDir','reverse');
colormap(ax, palette); clim(ax, limits);
title(ax, name, 'Interpreter','none'); xlabel(ax, 'x pixel'); ylabel(ax, 'y pixel');
local_style(ax, cfg, 'dataset');
end

function count = local_plot_targets(ax, C, E, cfg, mode, ra0, labels, halo)
hold(ax, 'on'); count = 0;
for j = 1:numel(E.Targets)
    T = E.Targets(j);
    [x, y, valid] = local_target_coordinates(T, mode, ra0);
    last = find(valid, 1, 'last'); if isempty(last), continue; end
    count = count+1; color = local_target_color(C, j);
    % Control the linewidth of the target
    if halo, plot(ax, x, y, '-', 'Color','w', 'LineWidth',1); end
    plot(ax, x, y, '-', 'Color',color, 'LineWidth',1);
    if strcmp(T.Appearance, 'object') || nnz(valid) == 1
        if halo, plot(ax, x(last), y(last), 'o', 'Color','w', 'MarkerSize',2, 'LineWidth',2.5); end
        plot(ax, x(last), y(last), 'o', 'Color',color, ...
            'MarkerFaceColor',color, 'MarkerSize',6, 'LineWidth',1.3);
    else
        % Control the marker size in the target
        first = find(valid, 1, 'first');
        plot(ax, x(first), y(first), '.', 'Color',color, 'MarkerSize',2);
    end
    if labels
        text(ax, x(last), y(last), sprintf(' T%d', T.ID), 'Color',color, ...
            'FontWeight','bold', 'FontSize',cfg.TargetLabelFontSize, ...
            'FontName',cfg.FontName, 'Interpreter','none', 'Clipping','on');
    end
end
end

function [x, y, valid] = local_target_coordinates(T, mode, ra0)
H = T.Samples;
switch mode
    case 'camera', x = H.Xi_deg; y = H.Eta_deg;
    case 'radec', x = local_wrap(H.RA_deg-ra0); y = H.DEC_deg;
    otherwise, x = H.x; y = H.y;
end
x = reshape(double(x), 1, []); y = reshape(double(y), 1, []);
valid = reshape(H.Visible, 1, []) & reshape(H.RadiometricRate, 1, []) > 0 & ...
    isfinite(x) & isfinite(y);
x(~valid) = NaN; y(~valid) = NaN;
% NaN gaps prevent drawing through occultation, eclipse or the RA wrap seam.
if strcmp(mode, 'radec')
    seam = [false abs(diff(x)) > 180];
    x(seam) = NaN; y(seam) = NaN; valid(seam) = false;
end
end

function color = local_target_color(C, j)
palette = C.Plot.Colors;
color = palette(1+mod(j-1,size(palette,1)),:);
end

function local_draw_boxes(ax, C, E, cfg)
hold(ax, 'on');
for j = 1:numel(E.Targets)
    T = E.Targets(j); if ~T.HasAnnotation, continue; end
    b = T.BoundingBox_xywh;
    if any(~isfinite(b)) || any(b(3:4) <= 0), continue; end
    color = local_target_color(C, j);
    % Stored xywh uses one-based pixel centres; draw edges at half pixels.
    rectangle(ax, 'Position',[b(1:2)-.5 b(3:4)], 'EdgeColor',color, 'LineWidth',1);
    text(ax, b(1), max(1,b(2)-3), sprintf('T%d',T.ID), 'Color',color, ...
        'FontWeight','bold', 'FontSize',max(cfg.DatasetLabelFontSize-1,4.5), ...
        'Interpreter','none', 'VerticalAlignment','bottom', ...
        'BackgroundColor','black', 'Margin',.5, 'Clipping','on');
end
end

function message = local_summary(C, E, cfg)
R = E.Targets;
lines = {sprintf('Exposure: %d (%.4g s)', E.Index, seconds(E.EndUTC-E.StartUTC)), ...
    ['Scene: ' local_case_tag(struct('CaseName',''), E)], ...
    sprintf('Image: %d x %d px', C.Camera.Nx, C.Camera.Ny), ...
    sprintf('Annotated targets: %d / %d', nnz([R.HasAnnotation]), numel(R)), ...
    sprintf('Object-mask pixels: %d', nnz(E.Labels.Semantic == 1)), ...
    sprintf('Trail-mask pixels: %d', nnz(E.Labels.Semantic == 2)), ...
    ' ', 'ID  type    bbox [xmin ymin xmax ymax]'};
rows = find([R.HasAnnotation]); n = min(numel(rows), cfg.DatasetMaxSummaryRows);
for j = 1:n
    T = R(rows(j)); b = T.BoundingBox_xywh;
    lines{end+1} = sprintf('%-3d %-7s [%d %d %d %d]', ...
        T.ID, T.Appearance, b(1), b(2), b(1)+b(3)-1, b(2)+b(4)-1); %#ok<AGROW>
end
if numel(rows) > n, lines{end+1} = sprintf('... %d additional targets', numel(rows)-n); end
message = strjoin(lines, newline);
end

function q = local_star_indices(count, limit)
countShown = min(count, limit);
if countShown == 0, q = []; else, q = unique(round(linspace(1,count,countShown))); end
end

function a = local_wrap(a)
a = mod(a+180, 360)-180;
end

function limits = local_limits(a, percentiles)
v = sort(double(a(isfinite(a))));
if isempty(v), limits = [0 1]; return; end
if numel(v) > 1
    q = 1+(numel(v)-1)*percentiles/100;
    limits = interp1(1:numel(v), v, q, 'linear');
else
    limits = [v(1) v(1)];
end
if limits(2) <= limits(1), limits = [v(1) v(end)]; end
if limits(2) <= limits(1), limits(2) = limits(1)+max(1,abs(limits(1)))*1e-6; end
end

function local_single_title(ax, cfg, E, label)
if cfg.ShowSingleAxisTitles
    title(ax, sprintf('%s%s | %s', cfg.TitlePrefix, label, local_case_tag(cfg,E)), ...
        'Interpreter','none');
end
end

function local_style(ax, cfg, mode)
switch mode
    case 'stage'
        axisSize = cfg.StageAxisFontSize; labelSize = cfg.StageLabelFontSize;
        titleSize = cfg.StagePanelTitleFontSize;
    case 'dataset'
        axisSize = cfg.DatasetAxisFontSize; labelSize = cfg.DatasetLabelFontSize;
        titleSize = cfg.DatasetPanelTitleFontSize;
    otherwise
        axisSize = cfg.AxisFontSize; labelSize = cfg.LabelFontSize;
        titleSize = cfg.LabelFontSize;
end
set(ax, 'Box','on', 'LineWidth',cfg.AxisLineWidth, 'FontSize',axisSize, ...
    'FontName',cfg.FontName, 'XColor','k', 'YColor','k');
ax.XLabel.FontSize = labelSize; ax.YLabel.FontSize = labelSize;
ax.Title.FontSize = titleSize;
if ~strcmp(mode, 'dataset'), grid(ax, 'on'); end
if cfg.ApplyTightInset && ~isa(ax.Parent, 'matlab.graphics.layout.TiledChartLayout')
    % Tiled layouts manage their own margins; single panels use TightInset.
    drawnow; set(ax, 'LooseInset',get(ax,'TightInset'));
end
end

function local_compact(ax, cfg)
if isfield(cfg, 'CompactNormalAxisLabels') && cfg.CompactNormalAxisLabels
    xlabel(ax, ''); ylabel(ax, '');
end
end

function local_panel_label(ax, label, cfg)
if ~cfg.AddPanelLabels, return; end
text(ax, .035, .965, label, 'Units','normalized', 'Color','w', ...
    'HorizontalAlignment','left', 'VerticalAlignment','top', ...
    'FontWeight','bold', 'FontSize',cfg.PanelLabelFontSize, 'Interpreter','none');
end
