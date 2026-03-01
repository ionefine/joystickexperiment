function [d, p, prms] = dynamic_contrast(d, p, prms)
% DYNAMIC_CONTRAST Fit the dynamic-contrast model to one dataset.
%
% Inputs
%   d     - analysis data struct (dc format)
%   p     - model parameter struct
%   prms  - fitting options struct
%
% Output
%   d, p, prms with fitted model outputs/parameters

% Per-sample validity mask (run x time)
d.ind = ones(size(d.response));

% Optionally keep only joystick-controlled trials.
if isfield(prms, 'joystickonly') && prms.joystickonly
    if ~isfield(d, 'joyused')
        warning('prms.joystickonly=true but d.joyused is missing. Keeping all runs.');
    else
        d.ind(~d.joyused, :) = 0;
    end
end

if isfield(prms, 'debuggingPlots') && prms.debuggingPlots
    % Show raw data (contrast + response) before cleaning.
    prms = dc.image_data(d, prms);
end

% Remove low-information runs.
d = dc.clean_data(d, prms, {'response_range'});

% Fit temporal kernel parameters (delay/tau) on binocular-only samples.
[d_bino, prms] = dc.select_data(d, 'bino', prms);
if isempty(d_bino.response)
    warning('No binocular samples found after cleaning; skipping hdr_delay/hdr_tau fit.');
else
    d_bino.prediction3D = d_bino.contrast;
    p = fitUW('dc.fit_model', p, {'hdr_delay', 'hdr_tau'}, prms, d_bino);
end

% Fit integration parameter(s) on all valid samples.
nFit = 1;
if isfield(prms, 'fit') && isnumeric(prms.fit) && isscalar(prms.fit) && prms.fit > 0
    nFit = round(prms.fit);
end

for i = 1:nFit %#ok<NASGU>
    d.prediction3D = d.contrast;
    p = fitUW('dc.fit_model', p, {'w'}, prms, d);
    [~, ~, d] = dc.fit_model(p, prms, d);
end
