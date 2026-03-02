% JOYSTICKEXPERIMENT_ANALYSIS
%
% Batch-fit the dynamic-contrast model for JoystickExperiment output files.

clear;
close all;

% Optional UW toolbox path (required for fitUW in many setups).
uwToolboxPath = 'C:\Users\Ione Fine\OneDrive - UW\Documents\code\toolboxes\UWToolbox\';
if isfolder(uwToolboxPath)
    addpath(genpath(uwToolboxPath));
end

% ------------------------- Analysis options -------------------------
dataType = 'psychophysics_bandpass';
stimulusType = 'congruent';
defaultOutputRoot = 'output';

% Let the user choose the output folder to analyse.
outputRoot = je.uigetdirOnTop(defaultOutputRoot, 'Select folder that contains participant output files');
if isequal(outputRoot, 0)
    error('Folder selection was cancelled.');
end

% Discover files created by JoystickExperiment_Run:
%   output/<subjectId>/<subjectId>_congruent_psychophysics_bandpass.mat
filePattern = sprintf('*_%s_%s.mat', stimulusType, dataType);
stimulusFiles = jea.find_stim_output_files(outputRoot, filePattern);

if isempty(stimulusFiles)
    error('No stimulus files found at %s using pattern %s', outputRoot, filePattern);
end

% Let user choose whether to analyse all files or a subset.
selectionMode = je.questdlgOnTop(sprintf('Found %d file(s). Analyse all or choose a subset?', numel(stimulusFiles)), ...
    'Select participant files', 'All files', 'Choose subset', 'Cancel', 'All files');

if isempty(selectionMode) || strcmp(selectionMode, 'Cancel')
    error('File selection was cancelled.');
elseif strcmp(selectionMode, 'Choose subset')
    fileLabels = strcat({stimulusFiles.folder}', filesep, {stimulusFiles.name}');
    [selectedFileIndex, ok] = je.listdlgOnTop('PromptString', 'Select participant files to analyse', ...
                                      'ListString', fileLabels, ...
                                      'SelectionMode', 'multiple', ...
                                      'ListSize', [900, 300]);
    if ~ok || isempty(selectedFileIndex)
        error('No files were selected for analysis.');
    end
    stimulusFiles = stimulusFiles(selectedFileIndex);
end

% Initialize result table lazily.
summaryTable = table();

% Timestamp used in output filename.
dateTimeString = datestr(now, 'yymmddHHMM');

for fileIndex = 1:numel(stimulusFiles)
    % Reset baseline params for each subject/file.
    [modelParams, fitOptions] = jea.params(dataType);
    close all;

    % ------------------------- Load file -------------------------
    stimulusFilePath = fullfile(stimulusFiles(fileIndex).folder, stimulusFiles(fileIndex).name);
    loadedData = load(stimulusFilePath);

    if ~isfield(loadedData, 'stim')
        warning('Skipping %s (no ''stim'' variable)', stimulusFilePath);
        continue;
    end

    % Subject identifier is stored with the session metadata.
    subjectId = 'unknown';
    if isfield(loadedData, 'session') && isfield(loadedData.session, 'subjectId')
        subjectId = loadedData.session.subjectId;
    end

    % Convert stimulus-format data -> analysis struct.
    analysisData = jea.stim_to_jea_data(loadedData.stim, subjectId);

    % Set identifiers for plotting/logging.
    modelParams.group = {'NA'};
    modelParams.id = {subjectId};

    jea.disp(['analysing subject ', subjectId, ' ', analysisData.datatype], fitOptions.verbalflag);

    % If fixation timestamps are present, remove peri-fixation samples.
    if isfield(analysisData, 'rfix_ts')
        analysisData = jea.fix_removal(analysisData);
    end

    % Optional contrast flip for bootstrap-like checks.
    if isfield(fitOptions, 'flipcontrast') && fitOptions.flipcontrast
        analysisData.contrast = 1 - analysisData.contrast;
    end

    % Set modality-appropriate HDR defaults.
    switch lower(analysisData.datatype)
        case {'psychophysics', 'vep_psychophysics', 'bold_psychophysics', 'psychophysics_bandpass'}
            modelParams.hdr_tau = 10 * analysisData.dt;
            modelParams.hdr_n = 1;
            modelParams.hdr_delay = 0.5;
        case {'vep', 'vep_allchan'}
            modelParams.hdr_tau = 50/1000;
            modelParams.hdr_n = 5;
            modelParams.hdr_delay = 0;
        case {'bold'}
            modelParams.hdr_tau = 1.5;
            modelParams.hdr_n = 3;
            modelParams.hdr_delay = 2;
        otherwise
            warning('Unrecognized datatype %s; using defaults from jea.params', analysisData.datatype);
    end

    % ------------------------- Fit model -------------------------
    [analysisData, modelParams, fitOptions] = jea.dynamic_contrast(analysisData, modelParams, fitOptions);
    [~, fitErrorNoCost, analysisData] = jea.fit_model(modelParams, fitOptions, analysisData);

    % Normalize attenuation vector and compute summaries.
    maxK = max(modelParams.k);
    if isfinite(maxK) && maxK > 0
        modelParams.k = modelParams.k ./ maxK;
    else
        warning('Skipping k normalization because max(k) is non-positive or non-finite.');
    end
    modelParams.min_k = min(modelParams.k);
    modelParams.err_noCost = fitErrorNoCost;

    if isfield(modelParams, 'U2')
        % Eyes organized as [RE LE].
        if modelParams.min_k == modelParams.k(1)
            modelParams.U_FEnAE = modelParams.U3 + 1;
            modelParams.U_AEnFE = modelParams.U2 + 1;
        else
            modelParams.U_FEnAE = modelParams.U2 + 1;
            modelParams.U_AEnFE = modelParams.U3 + 1;
        end
    end

    % Mean traces and quick quality metric.
    analysisData = jea.get_mean_tc(analysisData);
    modelParams.var = var([analysisData.mean(1).response analysisData.mean(2).response]);

    [~, fitErrorNoCost, ~] = jea.fit_model(modelParams, fitOptions, analysisData);
    modelParams.corr = 1 - fitErrorNoCost;

    % Save one row per subject/file.
    if isfield(modelParams, 'w')
        modelParams.w = abs(modelParams.w);
    end
    if isempty(summaryTable)
        summaryTable = struct2table(modelParams);
    else
        summaryTable(end+1, :) = struct2table(modelParams);
    end

    if isfield(fitOptions, 'pauseon') && fitOptions.pauseon
        pause;
    end
end

% ------------------------- Save summary table -------------------------
if ~isempty(summaryTable)
    modelName = 'model';
    if isfield(fitOptions, 'modelname') && ~isempty(fitOptions.modelname)
        modelName = fitOptions.modelname;
    end
    outputFile = fullfile('data', sprintf('%s_%s_%s.xlsx', dataType, modelName, dateTimeString));
    outputDir = fileparts(outputFile);
    if ~isempty(outputDir) && ~isfolder(outputDir)
        mkdir(outputDir);
    end
    writetable(summaryTable, outputFile);
    fprintf('Wrote %s\n', outputFile);
else
    warning('No analyzable files found; no output table was written.');
end
