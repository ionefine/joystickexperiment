function converted = convert_old_data_format(inputPath, outputPath)
%CONVERT_OLD_DATA_FORMAT Convert legacy joystick data MAT files to new stimulus format.
%
% converted = convert_old_data_format(inputPath)
% converted = convert_old_data_format(inputPath, outputPath)
%
% Supports two usage modes:
%   1) Single file conversion:
%      inputPath  = path/to/file.mat
%      outputPath = optional output .mat path (default: overwrite input file)
%
%   2) Directory conversion (recursive):
%      inputPath  = path/to/folder
%      outputPath = optional output folder (default: same folder tree)
%
% The function looks for legacy fields and rewrites them into:
%   stim.data.contrast, stim.data.response, stim.data.t,
%   stim.data.Bino_ON, stim.data.Mono_ON, stim.temporal.fastEye
%
% Returns:
%   converted - cell array of output file paths that were converted.

if nargin < 2
    outputPath = '';
end

if isfolder(inputPath)
    converted = convert_folder(inputPath, outputPath);
else
    converted = {convert_file(inputPath, outputPath)};
end
end

function converted = convert_folder(inputFolder, outputFolder)
files = dir(fullfile(inputFolder, '**', '*.mat'));
converted = {};

for i = 1:numel(files)
    inFile = fullfile(files(i).folder, files(i).name);

    if isempty(outputFolder)
        outFile = inFile;
    else
        relDir = erase(files(i).folder, [inputFolder filesep]);
        targetDir = fullfile(outputFolder, relDir);
        if ~isfolder(targetDir)
            mkdir(targetDir);
        end
        outFile = fullfile(targetDir, files(i).name);
    end

    try
        converted{end+1} = convert_file(inFile, outFile); %#ok<AGROW>
    catch ME
        warning('Skipping %s: %s', inFile, ME.message);
    end
end
end

function outFile = convert_file(inFile, outFile)
if nargin < 2 || isempty(outFile)
    outFile = inFile;
end

S = load(inFile);

if isfield(S, 'stim') && isfield(S.stim, 'data') && ...
        isfield(S.stim.data, 'contrast') && isfield(S.stim.data, 'response')
    % Already in new format; preserve and return.
    if ~strcmp(outFile, inFile)
        save(outFile, '-struct', 'S');
    end
    return;
end

[stim, session] = build_new_stim_struct(S, inFile);

S.stim = stim;
if ~isfield(S, 'session')
    S.session = session;
end

save(outFile, '-struct', 'S');
end

function [stim, session] = build_new_stim_struct(S, inFile)
subjectId = infer_subject_id(S, inFile);

stim = struct();
stim.data = struct();
stim.temporal = struct();
session = struct('subjectId', subjectId);

source = [];
if isfield(S, 'd') && isstruct(S.d)
    source = S.d;
elseif isfield(S, 'stim') && isstruct(S.stim)
    source = S.stim;
else
    source = S;
end

stim.data.contrast = pick_field(source, {'contrast', 'Contrast'});
stim.data.response = pick_field(source, {'response', 'Response'});

if isempty(stim.data.contrast) || isempty(stim.data.response)
    error('Could not find legacy contrast/response fields in this file.');
end

stim.data.t = pick_field(source, {'t', 'time', 'Time'});
if isempty(stim.data.t)
    stim.data.t = 1:size(stim.data.response, 2);
end

stim.data.Bino_ON = pick_field(source, {'Bino_ON', 'bino_on', 'binoOn'});
if isempty(stim.data.Bino_ON)
    stim.data.Bino_ON = zeros(size(stim.data.response));
end

stim.data.Mono_ON = pick_field(source, {'Mono_ON', 'mono_on', 'monoOn'});
if isempty(stim.data.Mono_ON)
    stim.data.Mono_ON = zeros(size(stim.data.response));
end

stim.temporal.fastEye = pick_field(source, {'fastEye', 'fast_eye'});
if isempty(stim.temporal.fastEye)
    stim.temporal.fastEye = zeros(size(stim.data.response, 1), 1);
end

% Keep metadata when available.
if isfield(S, 'opts')
    stim.opts = S.opts;
end
if isfield(S, 'display')
    stim.display = S.display;
end
end

function val = pick_field(s, names)
val = [];
for i = 1:numel(names)
    if isfield(s, names{i})
        val = s.(names{i});
        return;
    end
end

% Some legacy files place fields under .data
if isfield(s, 'data') && isstruct(s.data)
    for i = 1:numel(names)
        if isfield(s.data, names{i})
            val = s.data.(names{i});
            return;
        end
    end
end
end

function subjectId = infer_subject_id(S, inFile)
if isfield(S, 'session') && isstruct(S.session) && isfield(S.session, 'subjectId')
    subjectId = S.session.subjectId;
    return;
end
if isfield(S, 'd') && isstruct(S.d) && isfield(S.d, 'subjectId')
    subjectId = S.d.subjectId;
    return;
end

[~, name] = fileparts(inFile);
tok = regexp(name, '^([^_]+)_', 'tokens', 'once');
if ~isempty(tok)
    subjectId = tok{1};
else
    subjectId = 'unknown';
end
end
