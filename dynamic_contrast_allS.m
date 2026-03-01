% DYNAMIC_CONTRAST_ALLS
%
% Batch-fit the dynamic contrast model for all stimulus output files.
% Refactor goals:
%   1) Use naming consistent with stimulus code (subjectId, stim.data.*, etc.)
%   2) Read the MAT files produced by JoystickExperiment_Run.m
%   3) Improve readability and add comments
%   4) Add basic bug checks/guards

clear;
close all;

% Optional UW toolbox path (required for fitUW in many setups).
uwToolboxPath = 'C:\Users\Ione Fine\OneDrive - UW\Documents\code\toolboxes\UWToolbox\';
if isfolder(uwToolboxPath)
    addpath(genpath(uwToolboxPath));
end

% ------------------------- Analysis options -------------------------
datatype = 'psychophysics_bandpass';
stimType = 'congruent';
outputRoot = 'output';

% Discover files created by JoystickExperiment_Run:
%   output/<subjectId>/<subjectId>_congruent_psychophysics_bandpass.mat
filePattern = sprintf('*_%s_%s.mat', stimType, datatype);
stimFiles = dc.find_stim_output_files(outputRoot, filePattern);

if isempty(stimFiles)
    error('No stimulus files found at %s using pattern %s', outputRoot, filePattern);
end

% Initialize a result table lazily.
P = table();

% Timestamp used in output filename.
datetimeStr = datestr(now, 'yymmddHHMM');

for iFile = 1:numel(stimFiles)
    % Reset baseline params for each subject/file.
    [p, prms] = dc.params(datatype);
    close all;

    % ------------------------- Load file -------------------------
    filePath = fullfile(stimFiles(iFile).folder, stimFiles(iFile).name);
    S = load(filePath);

    if ~isfield(S, 'stim')
        warning('Skipping %s (no ''stim'' variable)', filePath);
        continue;
    end

    % Prefer saved session.subjectId when available.
    subjectId = 'unknown';
    if isfield(S, 'session') && isfield(S.session, 'subjectId')
        subjectId = S.session.subjectId;
    else
        tok = regexp(stimFiles(iFile).name, '^([^_]+)_', 'tokens', 'once');
        if ~isempty(tok)
            subjectId = tok{1};
        end
    end

    % Convert stimulus-format data -> dc-format analysis struct.
    d = dc.stim_to_dc_data(S.stim, subjectId);

    % Add identifiers expected by legacy downstream code.
    p.group = {'NA'};
    p.id = {subjectId};

    dc.disp(['analysing subject ', subjectId, ' ', d.datatype], prms.verbalflag);

    % If fixation timestamps are present, remove peri-fixation samples.
    if isfield(d, 'rfix_ts')
        d = dc.fix_removal(d);
    end

    % Optional contrast flip for bootstrap-like checks.
    if prms.flipcontrast
        d.contrast = 1 - d.contrast;
    end

    % Set modality-appropriate HDR defaults.
    switch lower(d.datatype)
        case {'psychophysics', 'vep_psychophysics', 'bold_psychophysics', 'psychophysics_bandpass'}
            p.hdr_tau = 10 * d.dt;
            p.hdr_n = 1;
            p.hdr_delay = 0.5;
        case {'vep', 'vep_allchan'}
            p.hdr_tau = 50/1000;
            p.hdr_n = 5;
            p.hdr_delay = 0;
        case {'bold'}
            p.hdr_tau = 1.5;
            p.hdr_n = 3;
            p.hdr_delay = 2;
        otherwise
            warning('Unrecognized datatype %s; using defaults from dc.params', d.datatype);
    end

    % ------------------------- Fit model -------------------------
    [d, p, prms] = dynamic_contrast(d, p, prms);
    [~, errNoCost, d] = dc.fit_model(p, prms, d);

    % Normalize attenuation vector and compute summaries.
    p.k = p.k ./ max(p.k);
    p.min_k = min(p.k);
    p.err_noCost = errNoCost;

    if isfield(p, 'U2')
        % Eyes organized as [RE LE] in this legacy interpretation.
        if p.min_k == p.k(1)
            p.U_FEnAE = p.U3 + 1;
            p.U_AEnFE = p.U2 + 1;
        else
            p.U_FEnAE = p.U2 + 1;
            p.U_AEnFE = p.U3 + 1;
        end
    end

    % Mean traces and quick quality metric.
    d = dc.get_mean_tc(d);
    p.var = var([d.mean(1).response d.mean(2).response]);

    [~, errNoCost, ~] = dc.fit_model(p, prms, d);
    p.corr = 1 - errNoCost;

    % Save one row per subject/file.
    p.w = abs(p.w);
    if isempty(P)
        P = struct2table(p);
    else
        P(end+1, :) = struct2table(p);
    end

    if prms.pauseon
        pause;
    end
end

% ------------------------- Save summary table -------------------------
if ~isempty(P)
    outputFile = fullfile('data', sprintf('%s_%s_%s.xlsx', datatype, prms.modelname, datetimeStr));
    outDir = fileparts(outputFile);
    if ~isempty(outDir) && ~isfolder(outDir)
        mkdir(outDir);
    end
    writetable(P, outputFile);
    fprintf('Wrote %s\n', outputFile);
else
    warning('No analyzable files found; no output table was written.');
end
