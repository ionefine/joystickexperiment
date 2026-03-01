function plot_params(mode, varargin)
%PLOT_PARAMS Unified plotting utility for model-fit parameters.
%   plot_params() plots all supported datasets (MRI, psychophysics, VEP,
%   and MRI ROI-pair/parameter relationships) from .xlsx tables in /data.
%
%   plot_params(MODE) runs one mode only:
%     'mri'            - grouped ROI summaries for MRI tables
%     'psychophysics'  - grouped summaries for psychophysics tables
%     'vep'            - grouped summaries for VEP tables
%     'mri_individual' - scatter/correlation views across MRI ROIs
%
%   Optional name/value pairs:
%     'DataDir'  (default: fullfile(pwd,'data'))
%     'GroupList' (default: {'NS','BD','AM'})
%     'Verbose' (default: true)
%
% Notes:
% - This file replaces the old split scripts:
%   plot_params_mri.m, plot_params_psycho.m, plot_params_vep.m,
%   plot_params_mri_ind.m.
% - Compatible with newer stimulus naming conventions by normalizing labels
%   like "congruent_psychophysics_bandpass" -> "congruent".

if nargin < 1 || isempty(mode)
    mode = 'all';
end

opts = parseInputs(varargin{:});
mode = lower(string(mode));

switch mode
    case "all"
        runGroupedMode('mri', opts);
        runGroupedMode('psychophysics', opts);
        runGroupedMode('vep', opts);
        runMriIndividualMode(opts);
    case {"mri","psychophysics","vep"}
        runGroupedMode(char(mode), opts);
    case "mri_individual"
        runMriIndividualMode(opts);
    otherwise
        error('Unknown mode "%s".', mode);
end
end

function opts = parseInputs(varargin)
p = inputParser;
p.addParameter('DataDir', fullfile(pwd, 'data'), @(x) ischar(x) || isstring(x));
p.addParameter('GroupList', {'NS','BD','AM'}, @(x) iscellstr(x) || isstring(x));
p.addParameter('Verbose', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
opts = p.Results;
opts.DataDir = char(opts.DataDir);
opts.GroupList = cellstr(opts.GroupList);
end

function runGroupedMode(modality, opts)
files = findTables(opts.DataDir, modality);
if isempty(files)
    if opts.Verbose
        warning('plot_params:%s:noFiles', modality, 'No %s tables found in %s.', modality, opts.DataDir);
    end
    return;
end

paramList = {'err_noCost','min_k','w'};
for i = 1:numel(files)
    T = readtable(files{i});
    T = normalizeTable(T);

    if ~all(ismember({'group','loc'}, T.Properties.VariableNames))
        warning('plot_params:%s:missingColumns', modality, 'Skipping %s (missing group/loc columns).', files{i});
        continue;
    end

    titleStub = sprintf('%s: %s', upper(modality), stripExtension(files{i}));

    switch modality
        case 'mri'
            if any(strcmp(T.Properties.VariableNames, 'loc'))
                xCats = unique(T.loc, 'stable');
            else
                xCats = {'V1','V2','V3'};
            end
            plotGrouped(T, paramList([1 3]), xCats, opts.GroupList, titleStub, true, true);
        case {'psychophysics','vep'}
            if any(strcmp(T.Properties.VariableNames, 'loc'))
                firstLoc = T.loc{1};
                T = T(strcmp(T.loc, firstLoc), :);
            end
            plotGrouped(T, paramList, opts.GroupList, opts.GroupList, titleStub, false, false);
    end
end
end

function runMriIndividualMode(opts)
files = findTables(opts.DataDir, 'mri');
if isempty(files)
    if opts.Verbose
        warning('plot_params:mri_individual:noFiles', 'No MRI tables found in %s.', opts.DataDir);
    end
    return;
end

T = normalizeTable(readtable(files{1}));
required = {'group','loc','min_k','w'};
if ~all(ismember(required, T.Properties.VariableNames))
    warning('plot_params:mri_individual:missingColumns', 'Missing required columns in %s.', files{1});
    return;
end

roiList = unique(T.loc, 'stable');
colors = [0.2 1 0.2; 0 0.7 0.1; 0 0.3 0];

% ROI-wise min_k vs w relation
figure('Name', 'MRI individual: min_k vs w'); clf;
tiledlayout(1, numel(roiList), 'Padding', 'compact', 'TileSpacing', 'compact');
for r = 1:numel(roiList)
    nexttile;
    R = T(strcmp(T.loc, roiList{r}), :);
    scatterByGroup(R, 'min_k', 'w', opts.GroupList, colors);
    lsline;
    axis square;
    xlim([0 1]); ylim([0 1]);
    xlabel('attenuation (k)'); ylabel('binocular integration (w)');
    title(roiList{r}, 'Interpreter', 'none');
end

% Cross-ROI consistency between first two ROIs (if available)
if numel(roiList) >= 2 && any(strcmp(T.Properties.VariableNames, 'id'))
    figure('Name', 'MRI individual: ROI pair consistency'); clf;
    tiledlayout(1,2, 'Padding', 'compact', 'TileSpacing', 'compact');
    for p = 1:2
        nexttile;
        hold on;
        for g = 1:numel(opts.GroupList)
            G = T(strcmp(T.group, opts.GroupList{g}), :);
            [xv, yv] = pairById(G, roiList{1}, roiList{2}, ternary(p==1, 'min_k', 'w'));
            scatter(xv, yv, 22, colors(g,:), 'filled', 'MarkerFaceAlpha', .55, 'MarkerEdgeAlpha', .55);
        end
        plot([0 1],[0 1],'k-');
        axis square;
        xlim([-0.05 1.05]); ylim([-0.05 1.05]);
        xlabel(roiList{1}); ylabel(roiList{2});
        title(ternary(p==1, 'k', 'w'));
    end
end
end

function plotGrouped(T, params, xCategories, groupList, titleStub, absValues, xIsLoc)
colors = [0.2 1 0.2; 0 0.7 0.1; 0 0.3 0];
lims = struct('err_noCost',[0 0.13], 'min_k',[-0.5 1.1], 'w',[-0.5 1.1]);

for p = 1:numel(params)
    param = params{p};
    if ~any(strcmp(T.Properties.VariableNames, param))
        continue;
    end

    figure('Name', sprintf('%s | %s', titleStub, param)); clf;
    hold on;

    gMask = true(height(T),1);
    if strcmp(param, 'err_noCost') && height(T) > 10
        gMask = ~isoutlier(T.(param));
    end

    for g = 1:numel(groupList)
        for x = 1:numel(xCategories)
            if xIsLoc
                xMask = strcmp(T.loc, xCategories{x});
            else
                xMask = strcmp(T.group, xCategories{x});
            end
            rows = strcmp(T.group, groupList{g}) & xMask & gMask;
            y = T.(param)(rows);
            y = y(~isoutlier(y));
            y = y(~isnan(y));
            if isempty(y)
                continue;
            end

            if absValues
                yPlot = abs(y);
            else
                yPlot = y;
            end

            xPos = x + (g-2)*0.08;
            scatter(xPos + 0.03*randn(size(yPlot)), yPlot, 16, colors(g,:), 'filled', ...
                'MarkerFaceAlpha', .5, 'MarkerEdgeAlpha', .5);
            med = median(yPlot, 'omitnan');
            q = prctile(yPlot, [25 75]);
            scatter(xPos, med, 44, 's', 'MarkerFaceColor', colors(g,:), ...
                'MarkerEdgeColor', colors(g,:), 'MarkerFaceAlpha', .6, 'MarkerEdgeAlpha', .6);
            line([xPos xPos], q, 'Color', colors(g,:), 'LineWidth', .75);
        end
    end

    xlim([0.5 numel(xCategories)+0.5]);
    set(gca, 'XTick', 1:numel(xCategories), 'XTickLabel', xCategories);
    if isfield(lims, param)
        ylim(lims.(param));
    end
    title(titleStub, 'Interpreter', 'none');
    ylabel(param, 'Interpreter', 'none');
    legend(groupList, 'Location', 'best');
    box on;
end
end

function scatterByGroup(T, xName, yName, groupList, colors)
hold on;
for g = 1:numel(groupList)
    G = T(strcmp(T.group, groupList{g}), :);
    if isempty(G), continue; end
    scatter(G.(xName), G.(yName), 24, colors(g,:), 'filled', ...
        'MarkerFaceAlpha', .55, 'MarkerEdgeAlpha', .55);
end
end

function [xVals, yVals] = pairById(T, locA, locB, param)
ids = unique(T.id, 'stable');
xVals = nan(numel(ids),1);
yVals = nan(numel(ids),1);
for i = 1:numel(ids)
    I = T(strcmp(T.id, ids{i}), :);
    a = I.(param)(strcmp(I.loc, locA));
    b = I.(param)(strcmp(I.loc, locB));
    if ~isempty(a) && ~isempty(b)
        xVals(i) = a(1);
        yVals(i) = b(1);
    end
end
good = ~(isnan(xVals) | isnan(yVals));
xVals = xVals(good);
yVals = yVals(good);
end

function T = normalizeTable(T)
% Normalize variable names and newer stimulus labels.
T.Properties.VariableNames = matlab.lang.makeUniqueStrings(lower(T.Properties.VariableNames));

% Canonical aliasing for common naming variants
T = alias(T, 'roi', 'loc');
T = alias(T, 'region', 'loc');
T = alias(T, 'groupname', 'group');

if any(strcmp(T.Properties.VariableNames, 'stim'))
    T.stim = cellfun(@normalizeStim, cellstr(string(T.stim)), 'UniformOutput', false);
end
end

function T = alias(T, fromName, toName)
if any(strcmp(T.Properties.VariableNames, fromName)) && ~any(strcmp(T.Properties.VariableNames, toName))
    T.(toName) = T.(fromName);
end
end

function out = normalizeStim(in)
out = lower(strtrim(in));
out = regexprep(out, '_psychophysics_bandpass$', '');
out = regexprep(out, '_bandpass$', '');
out = regexprep(out, '_new$', '');
end

function files = findTables(dataDir, modality)
if ~isfolder(dataDir)
    files = {};
    return;
end
allFiles = dir(fullfile(dataDir, '*.xlsx'));
names = lower(string({allFiles.name}));
switch modality
    case 'mri'
        keep = contains(names, "bold") & ~contains(names, "psycho") & ~contains(names, "vep");
    case 'psychophysics'
        keep = contains(names, "psycho");
    case 'vep'
        keep = contains(names, "vep");
    otherwise
        keep = false(size(names));
end
files = arrayfun(@(d) fullfile(d.folder, d.name), allFiles(keep), 'UniformOutput', false);
end

function s = stripExtension(pathStr)
[~, s, ~] = fileparts(pathStr);
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
