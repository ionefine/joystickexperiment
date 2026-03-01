function [d,p,prms] = dynamic_contrast(d, p, prms )%% full data set - initial data cleaning

d.ind = ones(size(d.response));

if prms.joystickonly     % 'true' then remove passive trials where joystick wasn't used
    d.ind(~d.joyused, :) = 0;
end

if prms.debuggingPlots
    prms = dc.image_data(d, prms); % s raw data: LE & RE contrast and response
end

d  = dc.clean_data(d, prms, {'response_range'}); % get rid of runs that didn't vary across a decent range during calibration

%% binocular dataset - find hdr, delay & response slope

%% monocular dataset - calculate monocular attenuation

[d_bino, prms] = dc.select_data(d, 'bino', prms);
d_bino.prediction3D = d_bino.contrast;
freeList = {'hdr_delay', 'hdr_tau'};
p = fitUW('dc.fit_model', p, freeList, prms, d_bino);


for i = 1:prms.fit

    % if sum(contains(prms.freeList, 'k')) &  strcmp(p.group, 'AM')
    %     %  old_errFtn = p.errFtn;
    %     %  p.errFtn = 'corr';
    %     [d_mono,prms] = dc.select_data(d, 'mono', prms);
    %     d_mono.prediction3D = d_mono.contrast;
    %     p = fitUW('dc.fit_model', p, {'k'}, prms, d_mono);
    %     p.k = p.k./max(p.k);
    %     % p.errFtn = old_errFtn;
    % end
    % no_k = ~contains(prms.freeList, 'k');
    %% fit & plot remainder of the parameters
    d.prediction3D = d.contrast;
    p = fitUW('dc.fit_model', p, {'w'}, prms, d);
    [errtn, errnoCost, d] = dc.fit_model(p, prms, d);
    %  p = fitUW('dc.fit_model', p, {'k'}, prms, d);
    % p.k = p.k./max(p.k);
    %  p = fitUW('dc.fit_model', p, {'w'}, prms, d);
    %  [errtn, errnoCost, d] = dc.fit_model(p, prms, d);


    % else % using preset parameters
    %     d.prediction3D = d.contrast;
    %     freeList = {'scalefac'};
    %     p = fitUW('dc.fit_model', p, freeList, prms, d);
    %     p.k = [1 1]; p.w = 0;
    %     [err, errnoCost, d] = dc.fit_model( p, prms, d);
    % end
end

