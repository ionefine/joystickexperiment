classdef dc
    methods(Static)
        function [d, Locs] = average_hemis(d, Locs, averagehemis)
            if averagehemis==true && (strcmp(d.datatype, 'vep_allchan') | strcmp(d.datatype, 'bold'))
                Locs =  Locs(Locs.average==true, :);
                for l = 1:size(Locs, 1)
                    for ll = 1:length(d.locs)
                        indx(ll, 1) = strcmp(d.locs(ll).labels, Locs.loc1(l));
                        indx(ll, 2) = strcmp(d.locs(ll).labels, Locs.loc2(l));
                    end
                    if length(find(indx)==2)
                        tmp(:, :, l) = mean(d.response(:, :, find(sum(double(indx), 2))),3);
                    else
                        tmp(:, :, l) = d.response(:, :, find(indx(:, 1)));
                    end
                end
                d = rmfield(d,'response');
                d.response = tmp;
            elseif strcmp(d.datatype, 'vep_allchan') | strcmp(d.datatype, 'bold')
                Locs = Locs(Locs.average==false, :);
                for l = 1:size(Locs, 1)
                    for ll = 1:length(d.locs)
                        indx(ll, 1) = strcmp(d.locs(ll).labels, Locs.loc1(l));
                        indx(ll, 2) = strcmp(d.locs(ll).labels, Locs.loc2(l));
                    end
                    tmp(:, :, l) = d.response(:, :, find(indx(:, 1)));
                end
                d = rmfield(d,'response');
                d.response = tmp;
            end
        end

        function [d, p, cost] = calc_model(prms, d, p)
            % Run through the  stages of the model
            cost = 0;
            for c = 1:length(prms.model_components)
                if strcmp(prms.model_components(c), 'monocular_attenuation')
                    [d, p, tmpCost] = dc.monocular_attenuation(prms, d, p);
                elseif strcmp(prms.model_components(c), 'conv_hdr')
                    [d, p, tmpCost] = dc.conv_hdr(prms, d, p);
                elseif strcmp(prms.model_components(c), 'normalization')
                    [d, p, tmpCost] = dc.normalization(prms, d, p);
                elseif strcmp(prms.model_components(c), 'weighted_meanmax')
                    [d, p, tmpCost] = dc.weighted_meanmax(prms, d, p);
                elseif strcmp(prms.model_components(c), 'weighted_binocular_norm')
                    [d, p, tmpCost] = dc.weighted_binocular_norm(prms, d, p);
                elseif strcmp(prms.model_components(c), 'weighted_monocular_norm')
                    [d, p, tmpCost] = dc.weighted_monocular_norm(prms, d, p);
                elseif strcmp(prms.model_components(c), 'binocular_norm')
                    [d, p, tmpCost] = dc.heeger_moradi(prms, d, p);
                elseif strcmp(prms.model_components(c), 'modifiedsoftmax')
                    [d, p, tmpCost] = dc.modifiedsoftmax(prms, d, p);
                elseif strcmp(prms.model_components(c), 'power_law')
                    [d, p, tmpCost] = dc.power_law(prms, d, p);
                elseif sum(ismember(prms.model_components(c), 'prctile_scale'))
                    [d, p, tmpCost] = dc.prctile_scale(prms, d, p);
                else
                    error('model component of prms.model_components not recognized');
                end
                cost = cost + tmpCost;
            end
        end

        function files = find_stim_output_files(rootDir, pattern)
            % Recursively find stimulus output MAT files under rootDir.
            %
            % Example:
            %   files = dc.find_stim_output_files('output', '*_congruent_psychophysics_bandpass.mat')
            if nargin < 2 || isempty(pattern)
                pattern = '*_congruent_psychophysics_bandpass.mat';
            end
            if nargin < 1 || isempty(rootDir)
                rootDir = 'output';
            end
            files = dir(fullfile(rootDir, '**', pattern));
        end

        function [d, subjectId] = stim_to_dc_data(stim, subjectId)
            % Convert saved stimulus output into the dc analysis format.
            %
            % Stimulus fields follow JoystickExperiment_Run / je.m naming:
            %   stim.data.contrast(run,time,eye)
            %   stim.data.response(run,time)
            %   stim.data.Bino_ON(run,time)
            %   stim.data.Mono_ON(run,time)
            %   stim.temporal.fastEye(run) [0=left fast, 1=right fast]

            if ~isfield(stim, 'data') || ~isfield(stim.data, 'contrast') || ~isfield(stim.data, 'response')
                error('Stim file missing required fields stim.data.contrast and/or stim.data.response.');
            end

            if nargin < 2 || isempty(subjectId)
                subjectId = 'unknown';
            end

            d = struct();
            d.datatype = 'psychophysics_bandpass';
            d.subjectId = subjectId;

            d.contrast = stim.data.contrast;
            d.response = stim.data.response;

            if isfield(stim.data, 't')
                d.t = stim.data.t;
                if numel(d.t) > 1
                    d.dt = d.t(2) - d.t(1);
                else
                    d.dt = 1;
                end
            else
                d.t = 1:size(d.response, 2);
                d.dt = 1;
            end

            if isfield(stim.data, 'Bino_ON')
                d.bino_on = double(stim.data.Bino_ON);
            else
                d.bino_on = zeros(size(d.response));
            end

            if isfield(stim.data, 'Mono_ON')
                d.mono_on = double(stim.data.Mono_ON);
            else
                d.mono_on = zeros(size(d.response));
            end

            if isfield(stim, 'temporal') && isfield(stim.temporal, 'fastEye')
                d.fastEye = stim.temporal.fastEye(:);
            else
                d.fastEye = zeros(size(d.response, 1), 1);
            end

            % New stimulus code always records joystick responses.
            d.joyused = true(size(d.response, 1), 1);
            d.goodruns = zeros(size(d.response, 1), 1);
        end
        function p = calc_response_slope(d, p)
            if size(d.prediction)~=size(d.response)
                error('must use meaned binocular data, contrast and response should be same sized vectors')
            else
                ind = find(~isnan(d.prediction.*d.response));
                p.linfit = polyfit(d.prediction(ind),d.response(ind), 1); % the linear scaling that best fits the data
            end
        end
        function d  = clean_data(d, prms, clean_style)
            % prms.cleanrange = (0-1); only calibrate or fit data where
            % there's this much range in the data (on a given trial)
            % d.response has 'clean' data only
            % n_range is a list of the runs kept
            %  d.goodruns = zeros(size(d.ind, 1),1);
            if contains(clean_style, 'response_range')
                %  disp('cleaning data based on inadequate response range');
                for r = 1:size(d.response, 1)
                    if max(d.response(r,:))-min(d.response(r, :))<=prms.clean_range
                        d.ind(r, :) = 0;
                        d.response(r,:) = NaN;
                    else
                        d.goodruns(r) = 1;  %  calibration used an adequate range
                    end
                end
                if prms.debuggingPlots
                    disp([num2str(sum(d.goodruns)), ' out of ', num2str(length(d.goodruns)), ' were kept'])
                end
            end
            if contains(clean_style, 'outliers_over_runs')
                tmp = 0;
                for e = 1:ndims(d.contrast)-1
                    for t = 1:size(d.contrast, 2)
                        if d.ind(1, t, 1)==1
                            tmp = tmp + var(d.contrast(:, t, e), [], 1);
                        end
                    end
                end
                if tmp>.001 % sometimes a tiny bit more than 0, tol error
                    error('trying to look for outliers over runs, but runs do not contain identical stimulus input');
                else
                    for t = 1:size(d.contrast, 2)
                        tf = isoutlier(d.response(:, t), "quartiles");
                        if sum(tf)>0
                            d.ind(find(tf), t) = 0;
                        end
                    end
                end
            end
            if contains(clean_style, 'noisesignal')
                figure(10); clf;
                tmp = diff(d.noiseresponse);
                hist(tmp(:)); drawnow; pause(2)
            end
        end
        function [d, p, cost] = conv_hdr(prms, d, p)
            % convolve with an HDR based on Boynton et al. '96
            if ~isfield(d, 'prediction2D') % early HDR before eyes combined, can either use single or eye-specific hdr values
                prediction = d.prediction3D;
            else
                prediction = d.prediction2D; % late HDR, after eyes combined, should use single hdr values
            end

            d = dc.gamma(d, p);
            if size(d.hdr, 2)>size(prediction, 3)
                errordlg('specifying eye-specific delay, but hdr is calculated after eyes combined');
            end
            prediction(isnan(prediction)) = 0.5;
            for e = 1:size(d.hdr, 2)
                for r = 1:size(prediction, 1)
                    tmp(r,:) = conv(prediction(r,:, e), d.hdr(:, e));
                end
                tmp = tmp(:, 1:size(prediction,2));
                tmp = max(prediction(:,:, e),[],'all').*tmp./max(tmp(:)); % rescale to original amplitude
                tmp(:, 1:length(d.hdr)) = NaN;
                if ~isfield(d, 'prediction2D')
                    d.prediction3D(:, :, e) = tmp;
                else
                    d.prediction2D = tmp;
                end
            end
            % zero out beginning, where the hdr hasn't reached asymptote
            d.contrast(:, 1:length(d.hdr), :) = NaN;
            d.response(1:length(d.hdr)) = NaN;
            cost =  dc.getCost(p.hdr_delay, prms.penalize_hdr_delay) + ...
                dc.getCost(p.hdr_tau, prms.penalize_hdr_tau) + dc.getCost_w(p.w) +dc.getCost(p.w, prms.penalize_w);
        end

        function disp(str, verbalflag)
            if verbalflag
                disp(str);
            end
        end

        function [err, err_noCost, d] = fit_model(p, prms, d)
            [d, p, cost] = dc.calc_model(prms, d, p);
            err_noCost = dc.getErr(d, p);
            err =  err_noCost + cost;
        end

        function [d, p] = fixed_nonlinearity(d, p)
            if ~isfield(p, 'm') || isnan(p.m)
                p.m = [1 1];
            end
            for i = 1:length(p.m)
                d.prediction2D(:, i) = d.prediction2D(:,i).^p.m(i);
            end
        end

        function d = fix_removal(d)
            % Nan's out data after a target event
            bw = 0.5; % the amount of time after keypress response to blank out
            for r = 1:size(d.response, 1)
                badvals = [d.t>d.rfix_ts(r,1) & d.t < d.rfix_ts(r,2)+bw];
                badvals = badvals + [d.t>d.rfix_ts(r,3) & d.t < d.rfix_ts(r,4)+bw];
                d.response(r, find(badvals)) = NaN;
            end
        end

        function [d, p, cost] = prctile_scale(prms, d, p)
            pc = prctile(nanmedian(d.response, 1), [5, 95]);
            d.prediction2D = dc.scale(d.prediction2D, pc(1), pc(2))+p.offset;
            cost = 0;
        end

        function cost = getCost(val, range)
            over = sum(100000 * abs(max(val - range(2), 0)).^4);
            under = sum(100000*abs(min(val - range(1),0)).^4);
            cost = sum(over(:) + under(:));
        end

        function cost = getCost_w(val)
            cost = max((1-val), 0).*.0001;
        end
        function d = gamma(d,p)
            %	y=Gamma(n,k,t)
            %	returns a gamma function on vector t
            %	y=(t/k).^(n-1).*exp(-t/k)/(k*factorial(n-1));
            %	which is the result of an n stage leaky integrator.

            %	6/27/95 Written by G.M. Boynton at Stanford University
            %   4/19/19 Simplified it for Psychology 448/538 at U.W.
            %
            for e = 1:length(p.hdr_delay) % 1 if same hdr for both eyes, 2 if eye-specific hdrs
                t = d.t - p.hdr_delay(e);
                y = (d.t/p.hdr_tau(e)).^(p.hdr_n-1).* ...
                    exp(-d.t/p.hdr_tau(e))/(p.hdr_tau(e)*factorial(p.hdr_n-1));
                y(t<0) = 0;
                hdr(:, e) = y;
            end

            cutoffCandidates = find(abs(diff(max(hdr,[], 2))) > .001);
            if isempty(cutoffCandidates)
                cutoff = size(hdr, 1);
            else
                cutoff = max(cutoffCandidates);
            end
            d.hdr=hdr(1:cutoff, :);
            peak = find(max(hdr,[], 2)==max(hdr(:)));
            d.hdrpeak = d.t(peak(1)); % the peak of the hdr, used to shift timecourses for visualization

        end
        function err = getErr(d, p)
            % Calculates MSE between data and model prediction

            if strcmp(p.errFtn, 'corr')
                for r = 1:1:size(d.response,1)
                    resp = d.response(r,:);
                    if strcmp(d.datatype, 'bold')
                        resp(1:6/d.dt) = NaN;
                    end
                    pred = d.prediction2D(r,:);
                    i = find(~isnan(pred.*resp));
                    if ~isempty(i)
                        tmp(r) = 1-corr(resp(i)',pred(i)');
                    else
                        tmp(r) = NaN;
                    end
                end
                err = nanmean(tmp);
            else
                err=0; ct = 0;
                for r = 1:1:size(d.response,1)
                    resp = d.response(r,:);
                    if strcmp(d.datatype, 'bold')
                        resp(1:6/d.dt) = NaN;
                    end
                    pred = d.prediction2D(r,:);

                    i = find(d.ind(r,:).*~isnan(resp).*~isnan(pred));
                    if strcmp(p.errFtn, 'rmse')
                        tmp = sqrt(sum((resp(i)-pred(i)).^2)); % add SSE to err
                    else
                        tmp = sum((resp(i)-pred(i)).^2); % add SSE to err
                    end

                    err = err + sum(tmp(:));
                    ct = ct + length(find(i));
                end
                if strcmp(p.errFtn, 'mse') || strcmp(p.errFtn, 'rmse')
                    err = err./ct;
                end
            end
        end
        function d = get_mean_tc(d)
            for e = 1:2 % for each 'task' (which eye faster
                fastEye = d.fastEye(1:size(d.mono_on, 1));
                mono_on = d.mono_on(fastEye+1==e,:);
                gvals = NaN(size(mono_on)); gvals(~mono_on) = 1;
                if contains(d.datatype, 'bold')
                    gvals = ones(size(gvals));
                end
                for c = 1:2
                    contrast = d.contrast(fastEye+1==e,:, c);
                    d.mean(e).contrast(:, c) = squeeze(nanmean(contrast.*gvals));
                end

                response = d.response(fastEye+1==e,:);
                d.mean(e).response = nanmean(response.*gvals);
                d.mean(e).response_ste = nanstd(response.*gvals)./sqrt(nansum(gvals, 1));
                if isfield(d, 'prediction2D')
                    prediction2D = d.prediction2D(fastEye+1==e,:);
                    d.mean(e).prediction2D = nanmean(prediction2D.*gvals);
                end
            end
        end
        function [pBest,errBest] = gridsearch(funName,p,gridParams,gridList, varargin)
            % [pBest,errBest] = gridsearch(funName,p,gridParams,gridList, var1, var2, ..)
            %
            % Grid search to find best initial parameters for optimization.
            %
            % Inputs:
            %   funName            name of error function. Must have form compatible
            %                      with 'fit' and 'fitcon' :
            %                       [err] = <funName>(params, var1, var2, ...)
            %   p                  Structure with initial parameters
            %   gridParams         List of parameters (fields of p) to be gridded
            %   gridList           cell array of 1-d grid values, in order
            %                      corresponding to gridParams
            %   var1, etc.         additional variables to be passed in to funName
            %
            % Outputs:
            %   pBest              parameter structure with lowest err in the grid
            %   errBest            corresponding error value
            %
            % Note, parameters in 'p' that are not in 'gridParams' are kept at these
            % initial values for all function evaluations.
            %
            % Note also: gridParams can contain elements of an array, e.g. 'a(1)'.

            nParams = length(gridParams);

            % Generate a string to be evaluated that has the form (depending on
            % nParams): '[X{1},X{2},X{3}] = ndgrid(gridList{1},gridList{2},gridList{3});'

            leftStr = '[';
            rightStr = ' = ndgrid(';
            for i=1:nParams
                leftStr = strcat(leftStr,sprintf('X{%d},',i));
                rightStr = strcat(rightStr,sprintf('gridList{%d},',i));
            end
            str = sprintf('%s] %s);',leftStr(1:end-1),rightStr(1:end-1));
            eval(str);

            errBest = 1e10;
            pBest = p;

            % Loop through elements of the grid, saving the current best fit
            for i=1:length(X{1}(:))
                for j=1:nParams
                    % generate a string to be evaluated of the form (for example):
                    % 'p.a(1) = X{1}(1);'
                    str = sprintf('p.%s = X{%d}(%d);',gridParams{j},j,i);
                    eval(str)
                end

                % evaluate the function with these parameters
                err = feval(funName, p,varargin{:});

                % save the fit if it's currently the best
                if err<errBest
                    pBest = p;
                    errBest = err;
                end
            end
        end
        function prms = image_data(d, prms)
            figure(prms.plotNum); clf;
            set(gcf, 'Name', ['datatype: ', d.datatype]);

            subplot(3, 1, 1)
            imagesc(d.contrast(:, : ,1));
            colormap(gray);
            title('LE contrast');
            ylabel('runs'); xlabel('time');

            subplot(3, 1, 2)
            imagesc(d.contrast(:, : ,2));
            colormap(gray);
            title('RE contrast');
            ylabel('runs'); xlabel('time');

            disp(['min contrast = ', num2str(min(d.contrast(:))), ' max contrast = ', num2str(max(d.contrast(:)))]);

            subplot(3, 1, 3)
            imagesc(d.response);
            colormap(gray);
            %yticks(1:size(d.response,1));
            disp(['min response = ', num2str(min(d.response(:))), ' max response = ', num2str(max(d.response(:)))]);
            title('response');
            ylabel('runs'); xlabel('time');
            prms.plotNum = prms.plotNum + 1;
        end
        function [kfoldErr_mean, kfoldErr_ste, kfoldErr, kf_p] = kfold(p,prms, d, freeList)
            if ~isfield(p, 'nKfolds')
                p.nKfolds = size(d.response,1)/2;
            end
            tt = randi(p.nKfolds, size(d.response,1), 1);
            for f = 1:p.nKfolds
                train_ind = find(tt~=f);
                test_ind = find(tt==f);

                % do the training
                kf_d = d;
                kf_d.response = d.response(train_ind, :);
                kf_d.contrast= d.contrast(train_ind, :, :);
                kf_d.prediction_c = NaN(size(kf_d.contrast));
                kf_p(f) = fit('dc.fit_model',p, freeList, prms, kf_d);

                % grab error for training set
                kf_d = d;
                kf_d.response = d.response(test_ind, :);
                kf_d.contrast= d.contrast(test_ind, :, :);
                kf_d.prediction_c = NaN(size(kf_d.contrast));
                [err, err_noCost, kf_d] = dc.fit_model(kf_p(f), prms, kf_d);
                kfoldErr(f) = err_noCost;
            end
            kfoldErr_mean = nanmean(kfoldErr);
            kfoldErr_ste = nanstd(kfoldErr)./sqrt(length(kfoldErr));
            if find(isnan(kfoldErr))
                warning('kfoldErr is producing Nans!')
            end
        end
        function [d, p] = linscale(d, p)
            p = dc.calc_response_slope(d, p);
            [d, p] = dc.scale_response_slope(d, p);
        end
        function logy2raw(base, precision)
            % logy2raw(base, precision)
            %
            % Converts Y-axis labels from log to raw values.
            %
            % Inputs:
            %   base           Base of log transform (default: e)
            %   precision      Number of decimal places (default: 2)
            %
            % Example:
            % x = linspace(-3,0,11);
            % plot(log(x), log(x.^2));
            % logx2raw();
            % logy2raw(); % should be tolerant to multiple calls
            %
            % Note:
            % - See also: logx2raw.m

            % 11/17/96       gmb wrote it.
            % 6/6/96	     gmb added precision argument
            % 01/30/02       gmb updated it to use cell arrays, and to use original
            %                xtick values instead of converting labels. This way,
            %                multiple calls to this function doesn't keep converting
            %                the axis.
            % Edited by Kelly Chang - February 18, 2017

            %% Input Control

            if ~exist('base', 'var')
                base = exp(1);
            end

            if ~exist('precision', 'var')
                precision = 2;
            end

            %% Calculate Log x-axis Labels

            precision = sprintf('%%%2.1ff', precision*1.1);
            origYTick = get(gca, 'YTick'); % log y-axis labels (raw)
            newYTick = base.^(origYTick); % convert to raw
            newYLabel = arrayfun(@(x) sprintf(precision,x), newYTick, ...
                'UniformOutput', false); % write new y-axis labels
            set(gca, 'YTickLabel', newYLabel); % set y-axis labels of current graph
        end
        function  d = mean_data(d)
            v = var(d.contrast, [], 1); % making sure averaging over identical data
            if sum(d.ind)<sum(d.ind.*double(v<0.001))
                error(['trying to average over runs containing different contrasts ', prms.eyelabel{e}]);
            else
                d.contrast = nanmean(d.contrast, 1);
                d.response = nanmean(d.response, 1);
                d.ind = logical(max(d.ind, [], 1));
                if isfield(d, 'prediction')
                    d.prediction = nanmean(d.prediction, 1);
                end
            end
        end
        function [d, p, cost] = monocular_attenuation(prms, d, p)
            % linearly attenuates the contrast to each eye based on
            % parameter p.k

            if ~isfield(p, 'k'), p.k = [1 1]; end
            p.k = p.k./max(p.k);
            for e = 1:length(p.k)
                d.prediction3D(:, :, e) = p.k(e).*d.prediction3D(:,:, e);
            end
            cost = 0;
            for i = 1:length(p.k)
                cost = cost + dc.getCost(p.k(i), prms.penalize_k);
            end
            p.k = p.k./max(p.k);
        end
        function [d, p, cost] = normalization(prms, d, p)
            if ~isfield(p, 'sigma'); p.sigma = .1; end
            tmp(:, :,1) = d.prediction3D(:, :,1)./(p.U1*d.prediction3D(:, :,1) + p.U2*d.prediction3D(:, :,2)+p.sigma);
            tmp(:,:, 2) = d.prediction3D(:, :,2)./(p.U3*d.prediction3D(:, :,1) + p.U4*d.prediction3D(:, :,2)+p.sigma);
            cost = 0;
            cost =  dc.getCost(p.U1, prms.penalize_U) + dc.getCost(p.U2, prms.penalize_U) + ...
                dc.getCost(p.U3, prms.penalize_U) +dc.getCost(p.U4, prms.penalize_U) ;
            d.prediction3D = tmp;
        end
        function [p, prms] = params(datatype)
            % flags
            prms.savePlotOn = 0;        % if 1, saves plots
            prms.plotNum = 1;
            prms.debuggingPlots = 0;         % puts up a series of debugging plots
            prms.pauseon = 0;
            prms.averagehemis = true;   % average bold and vep data across hemispheres
            prms.flipcontrast = 0;
            prms.verbalflag = 0;

            p.k = [1 1];
            p.scalefac = 1; p.offset = 0;
            prms.fit = 1;
            prms.modelname = 'weighted_binocular_norm';%'weighted_meanmax'; % %%; %';%'; %

            if strcmp(prms.modelname, 'weighted_binocular_norm')
                prms.model_components = {'monocular_attenuation',   'weighted_binocular_norm','conv_hdr','prctile_scale'};
                prms.freeList = { 'w'};   p.errFtn = 'mse';
                p.n1 = 1;  p.sigma1 = 10; p.n2 = 10;p.sigma2 = 15;  p.pl = 1;
                p.w = 1;
            elseif strcmp(prms.modelname, 'binocular_norm')
                prms.model_components = {'monocular_attenuation',   'binocular_norm',  'conv_hdr', 'prctile_scale'};
                prms.freeList = {'k', 'n', 'sigma'};
                p.sigma = .1; p.n = 3.5;
            elseif strcmp(prms.modelname, 'weighted_meanmax')
                prms.model_components = {'monocular_attenuation',   'weighted_meanmax',  'conv_hdr', 'prctile_scale'};
                prms.freeList = {'w', 'k'}; p.errFtn = 'mse';
                p.w = 1;
            elseif strcmp(prms.modelname, 'weighted_meannorm')
                prms.model_components = {'monocular_attenuation',   'weighted_meannorm',  'conv_hdr', 'prctile_scale'};
                prms.freeList = {'k', 'w'};
                p.w = .5;
            end
            if ~contains(datatype, 'bold')  && ~contains(datatype, 'psychophysics')% vep non-psychophysics data
                ind = contains(prms.model_components, 'conv_hdr');
                prms.model_components = prms.model_components(~ind);
            end

            %% eeg only
            prms.topographyParams = {'delay', 'w'}; % which parameters do you want to visualize on topography maps (eeg only)
            prms.joystickonly = 0;              % 0 all trials, 1 only keep trials where joystick used

            %%%%% prms details: %%%%%%
            prms.eyelabel = {'left', 'right'};
            prms.clean_range = .8; %  use runs with this min range in the data, cleanrange = 0, uses all data

            prms.penalize_w = [0 1]; %penalize weighted mean max values (w) outside this range
            prms.penalize_delay = [-1 1]; % penalize delays (p.delay) outside this range(sec)
            prms.penalize_k = [0 Inf]; % penalize monocular attenuation (p.k) outside this range
            prms.penalize_slope = [ 0 2]; % penalize calibration slopes outside this range
            prms.penalize_U = [ 0 Inf]; % penalize U values slopes outside this range
            prms.penalize_scalefac = [0 3];
            prms.penalize_sigma = [.001 .8];
            prms.penalize_pl = [.1 2];
            prms.penalize_n = [1 6];

            switch lower(datatype)
                case {'psychophysics' , 'vep_psychophysics' , 'bold_psychophysics','psychophysics_bandpass'}
                    prms.penalize_hdr_tau = [.4 2];
                    prms.penalize_hdr_delay = [.2 3];
                case {'vep'}
                    prms.penalize_hdr_tau = [.10/1000 .4/1000];
                    prms.penalize_hdr_delay = [0 .5];
                case {'bold'}
                    prms.penalize_hdr_tau = [ 1 2];
                    prms.penalize_hdr_delay = [ 2 5];
                case { 'vep_allchan'}
                    prms.penalize_hdr_tau = [.4 2];
                    prms.penalize_hdr_delay = [.2 3];
            end
            if prms.flipcontrast % don't fit k for flipped contrast
                ind = contains(prms.freeList, 'k');
                prms.freeList = prms.freeList(~ind);
            end
            prms.gridsearch.use = 0; % not sure this is working BUGCHECK
            prms.gridsearch.freeList = {'U(2)','U(3)','sigma'};
            prms.gridsearch.matrices = {0:0.25:3, 0:0.25:3, 0:0.1:1};   % list of grid vectors for each parameter in gridParams

        end
        function prms = plot_timecourses(d, prms, p)
            figure(prms.plotNum);
            clist = [1 0 0 .5; 0 0 1 .5];
            for e = 1:2 % for each eye
                subplot(2, 2, e)
                for c = 1:2 % for each contrast
                    plot(d.t, d.mean(e).contrast(:, c),'Color',clist(c,:), 'LineWidth',2); hold on
                end
                axis off
                xlim([0 max(d.t)]); ylim([-0.1 1.1]);
                if isfield(d, 'hdrpeak')
                    offset = round(mean([d.hdrpeak/d.dt length(d.hdr)]));
                else
                    offset = 1;
                end
                % if isfield(d.mean, 'prediction2D')
                %     plot(d.t(1:end-offset+1), d.mean(e).prediction2D(offset:end), 'g'); hold on
                % end
                subplot(2,2, e+2)
                if isfield(d.mean, 'prediction2D')
                    plot(d.t(1:end-offset+1), d.mean(e).prediction2D(offset:end), 'g'); hold on
                end
                shadedErrorBar(d.t(1:end-offset+1),d.mean(e).response(offset:end),d.mean(e).response_ste(offset:end), 'lineProps',{'k-','markerfacecolor','k'})
                xlim([0 max(d.t)]); ylim([-0.1 1.1]);
            end
            set(gcf, 'Name', [d.datatype, ' ',p.group{1}, ' ', p.id{1}]);
            prms.plotNum = prms.plotNum + 1;
        end
        function  prms = plot_contrast_response(d, prms)
            figure(prms.plotNum);  clf;
            set(gcf, 'Name', 'contrast response');
            for e = 1:2
                subplot(2, 1, 1)
                plot(d.mean(e).contrast(:, 1), d.mean(e).response, '.', 'Color', [e/3 e/3 e/3]); hold on
                plot([-.5 1.5], [-.5 1.5], 'k--');

                subplot(2, 1, 2)
                plot(d.contrast(:, :, 1), d.response, '.', 'Color', [e/3 e/3 e/3]); hold on;
                plot([-.5 1.5], [-.5 1.5], 'k--');
            end
            xlabel('contrast/prediction'); ylabel('response')
            set(gca, 'XLim', [-.5 1.5]);set(gca, 'YLim', [-.5 1.5]); prms.plotNum = prms.plotNum + 1;
        end
        function out = scale(im, varargin)
            %  function sim = scaleif(im, newMin, newMax, oldMin, oldMax)
            %  Scales an image such that its lowest value attains newMin and
            %  it’s highest value attains newMax.  OldMin and oldMax are not
            %  necessary but are useful when you don't want to use the true
            %  min or max value.
            %
            % Adapted by Ione Fine, based on code from Rick Anthony
            % 6/5/2000

            if nargin==1
                newMin=0; newMax=1; oldMin = min(im(:)); oldMax = max(im(:));
            elseif nargin==3
                newMin=varargin{1};newMax=varargin{2}; oldMin = min(im(:)); oldMax = max(im(:));
            elseif nargin==5
                newMin=varargin{1}; newMax=varargin{2}; oldMin=varargin{3}; oldMax=varargin{4};
            else
                disp('Wrong number of input arguments'); out = NaN;return
            end

            if newMin>=newMax
                disp('Sorry new min must be smaller than new max'); out = NaN;  return
            end
            if oldMin>=oldMax
                disp('Sorry old min must be smaller than old max');  out = NaN; return
            end
            delta = (newMax-newMin)/(oldMax-oldMin);
            out = delta*(im-oldMin) + newMin;
        end
        % function [d, p, cost] = scale_response(prms, d, p)
        %     d.prediction2D = p.offset+(1-p.scalefac)/2 + (d.prediction2D * p.scalefac);
        %     cost = dc.getCost(p.scalefac, prms.penalize_scalefac);
        % end
        function [d, prms] = select_data(d, str, prms)
            if strcmp(str, 'bino')
                ind2D = logical(d.bino_on.*d.ind);
            elseif strcmp(str, 'mono')
                ind2D = logical(d.mono_on.*d.ind);
            elseif strcmp(str, 'dich')
                ind2D = logical(~(d.mono_on+d.bino_on).*d.ind);
            end
            ind3D = cat(3, ind2D, ind2D);
            d.contrast(~ind3D) = NaN;
            if prms.debuggingPlots
                figure(prms.plotNum); clf;
                tmp = 2+(254*(d.contrast));
                cmap=[0 0 1; gray(255)];

                for e = 1:2
                    subplot(2, 1, e)
                    image(tmp(:, :, e)); colormap(cmap);
                    xlabel('time'); set(gcf, 'Name', str);
                end
                prms.plotNum = prms.plotNum + 1;
            end
        end

        % function d = scale_response_slope(d, p)
        %     d.prediction_r = (d.response-p.linfit(2))/p.linfit(1);
        % end
       
        function [d, p, cost] = monocular_norm(prms, d, p)
            % monocular normalization
            % each eye undergoes contrast normalization separately
            d.prediction3D(:, :, 1) = (d.prediction3D(:, :, 1).^p.n)./(d.prediction3D(:, :,  1).^p.n+p.sigma);
            d.prediction3D(:, :, 2) = (d.prediction3D(:, :, 2).^p.n)./(d.prediction3D(:, :,  2).^p.n+p.sigma);
            cost = dc.getCost(p.sigma, prms.penalize_sigma);
        end
        function [d, p, cost] = binocular_norm(prms, d, p)
            % binocular normalization (heeger style)
            d.prediction2D =  (d.prediction3D(:, :, 1).^p.n + d.prediction3D(:, :, 2).^p.n) ./ ...
                (d.prediction3D(:, :, 1).^p.n + d.prediction3D(:, :, 2).^p.n + p.sigma);

            cost = dc.getCost(p.sigma, prms.penalize_sigma);
        end
        function [d, p, cost] = weighted_meanmax(prms, d, p)
            % this is the "mixture" model, relative weighting of the mean
            % and the max rule
            d.prediction2D = (1-abs(p.w)) * (nanmean(d.prediction3D, 3))  + abs(p.w) * max(d.prediction3D, [], 3);
            cost = dc.getCost(p.w, prms.penalize_w);
            p.w = abs(p.w);    p.w = min(p.w, 1);
        end
        function [d, p, cost] = weighted_binocular_norm(prms, d, p)
            % this is the "mixture" model, relative weighting of the
            % monocular normalization and binocular normalization (Heeger
            % style)
            % tmp2_mean = nanmean(d.prediction3D, 3);
            % first do monocular contrast normalization

            % average the monocular response
            tmp1 = max(nanmean(d.prediction3D, 3));meanref = max(tmp1(:));
            tmp2 = max(d.prediction3D, [], 3);maxref = max(tmp2(:));

            d.prediction3D_M(:, :, 1) = p.sigma1*((d.prediction3D(:, :, 1).^p.n1)./(d.prediction3D(:, :,  1).^p.n1+p.sigma1));
            d.prediction3D_M(:, :, 2) = p.sigma1*((d.prediction3D(:, :, 2).^p.n1)./(d.prediction3D(:, :,  2).^p.n1+p.sigma1));
            tmp_mean = nanmean(d.prediction3D_M, 3);
       %     tmp_mean = dc.scale(tmp_mean, 0, meanref);

            % second stage of binocular normalization
            tmp_norm = p.sigma2* (d.prediction3D(:, :, 1).^p.n2 + d.prediction3D(:, :, 2).^p.n2) ./ ...
                (d.prediction3D(:, :, 1).^p.n2 + d.prediction3D(:, :, 2).^p.n2 + p.sigma2);
            tmp_norm = tmp_norm.^(1/p.n2);
     %       tmp_norm = dc.scale(tmp_norm, 0, maxref);

            d.prediction2D = ((1-abs(p.w)).*tmp_mean) + (abs(p.w)*tmp_norm);
            cost = dc.getCost(p.w, prms.penalize_w);
         %   cost = cost + .01*((p.w.*(1-min(p.k))).^8);
        end

        function [d, p, cost] = modifiedsoftmax(prms, d, p)
            % this is a softmax rule but goes between mean and max
            p = 1./(1-p.w);
            d.prediction2D = ((1-p.w+1) * sum(d.prediction3D.^p.w, 3)).^1/p.w;
            cost = dc.getCost(abs(p.w), prms.penalize_w);
            p.w = abs(p.w);
        end

        function [d, p, cost] = power_law(prms, d, p)
            % compressive nonlinearity
            d.prediction2D = d.prediction2D.^p.pl;
            cost = dc.getCost(p.pl, prms.penalize_pl);
            p.w = abs(p.w);
        end
    end
end
