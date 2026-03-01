%function dynamic_contrast()
%
% fits dynamic contrast data. The data is first pre-processed by Mark and
% then that data is passed through 'convert_data' so that eeg, psycho and
% fMRI all use metameric parameters.
clear
close

% curDir = pwd;
% cd('data');
% allfilesP = dir(['*', 'psychophysics', '*','.xlsx']);
% allfilesB = dir(['*', 'bold_weighted', '*','.xlsx']);
% 
% for ps = 1:length(allfilesP)
%    delete(allfilesP(ps).name);
%     for bd = 1:length(allfilesB)
%        delete(allfilesB(bd).name);
%     end
% end
% cd(curDir);

addpath(genpath('C:\Users\Ione Fine\OneDrive - UW\Documents\code\toolboxes\UWToolbox\'));
datatype_list = {'bold'}; %'bold_psychophysics', 'bold'}
%'bold'}; %,,  }; %{'psychophysics_bandpass'}; %, 'bold'}; %,  %{'vep_psychophysics'}; %; %'vep_kids'}; %bold'}; %'bold' 'vep','bold'}; %, vep_psychophysics', 'bold_psychophysics'}; %,'bold'};{}; % all possible datatypes
%% 
datetimestr = datestr(now, 'yymmddHHMM');
for dd = 1:length(datatype_list)
    clear P p
    T = readtable('sID_table.xlsx');
    T = T(strcmp(T.datatype, datatype_list{dd}), :); % select the datatype
    T = T(strcmp(T.stimtype, 'congruent'), :);
    %T = T(strcmp(T.group, 'AM'), :); % to limit to a particular group
    % T = T(strcmp(T.id, 'BA_15'), :); % to limit to a particular individual
    
    % load the locations
    [p, prms] = dc.params(datatype_list{dd});

    L = readtable('dc_location_table');
    L = L([strcmp(L.datatype, datatype_list{dd})], :);

    for i = 1:size(T, 1)
        [p, prms] = dc.params(datatype_list{dd});
        close all
        I = T(i,: );
        load(['data' filesep I.datatype{1}, filesep I.group{1} '_', I.id{1}, '_', I.stimtype{1}, '_', I.datatype{1}]);
        p.group = I.group; p.id = I.id; d.datatype = I.datatype{1};

        switch p.id{1}
            case 'BA_15' % aniso + strabn
                p.k = [1 .988];
            case 'IA_24'% aniso
                p.k = [1 .856]; 
            case 'ID_16' %strab
                p.k = [1 .470];
            case 'YA_26'  % aniso
                p.k = [.91 1];
            case 'CI_20'
                p.k = [1 1];
            case 'KR_19'
                p.k = [1 .8];
            case 'AA_17'
                p.k = [.979 1];
            case 'CW_25'
                p.k = [.833 1];
            case 'JX_19'
                p.k = [1 1];
            case 'PB_23'
                p.k = [1 .83];
            case 'PP_17'
                p.k = [1 .8];
            case 'XX_02'
                p.k = [1 .943];
            case 'XX_04'
                p.k = [.940 1];
            case 'XX_09'
                p.k = [1 1];
            case 'XX_10'
                p.k = [.962 1];
            case 'YN_17'
                p.k = [1 .993];
        end

        dc.disp(['analysing subject ', p.group{1},' ',  p.id{1}, ' ', d.datatype], prms.verbalflag);

        [d, Locs] = dc.average_hemis(d, L, prms.averagehemis);

        if isfield(d, 'rfix_ts') % if taking responses
            d = dc.fix_removal(d);
        end

        d_raw = d; % save the full set of responses
        % collect location labels & positions

        if prms.flipcontrast % flip the contrast in the  the two eyes, poor man's bootstrap
            d.contrast = 1-d.contrast;
        end

        for l = 1; %:size(Locs, 1)
            p.loc = Locs.label(l);
            p.X = Locs.x(l);
            p.Y = Locs.y(l);

            dc.disp(['fitting location ', p.loc{1}, ' / ', num2str(l), ' out of ', num2str(size(Locs, 1))], prms.verbalflag);
            %% reset fitted parameter values

            switch lower(d.datatype)
                case {'psychophysics' , 'vep_psychophysics' , 'bold_psychophysics','psychophysics_bandpass'}
                    p.hdr_tau = 10*d.dt;  p.hdr_n = 1; p.hdr_delay = 0.5; 
                case {'vep', 'vep_allchan'}
                    p.hdr_tau = 50/1000;  p.hdr_n = 5; p.hdr_delay = 0;
                case {'bold'}
                    p.hdr_tau = 1.5;  p.hdr_n = 3; p.hdr_delay = 2;
                case 'vep_kids'
                    p.age = I.age; p.hdr_tau = 50/1000;  p.hdr_n = 5; p.hdr_delay = 0;
            end

            % fit data
            d.response = d_raw.response(:, :, l);
            [d, p, prms] = dynamic_contrast(d, p, prms);
         
            [err, errnoCost, d] = dc.fit_model( p, prms, d);

            dc.disp(['w = ', num2str(p.w)], prms.verbalflag)

            p.k = p.k./max(p.k);
            p.min_k = min(p.k);
            p.err_noCost = errnoCost;

            dc.disp(['k = ', num2str(p.k(1)), ' ', num2str(p.k(2))], prms.verbalflag)
            dc.disp(['offset = ', num2str(p.offset)], prms.verbalflag)
            dc.disp(['error = ', num2str(p.err_noCost)], prms.verbalflag)
            dc.disp(['delay = ', num2str(p.hdr_delay)], prms.verbalflag)
            dc.disp(['tau = ', num2str(p.hdr_tau)], prms.verbalflag)

            if isfield(p, 'U2')
                % Eyes are organized as RE, LE
                if p.min_k == p.k(1) % RE amblyopic
                    p.U_FEnAE = p.U3 + 1; % normalization by the amblypic eye onto the fellow
                    p.U_AEnFE = p.U2 + 1;% normalization by the fellow eye onto the amblyopic
                else
                    p.U_FEnAE = p.U2 + 1;
                    p.U_AEnFE = p.U3 + 1;
                end
            end

            % plot and save plot
            d = dc.get_mean_tc(d);
            if  prms.debuggingPlots
                dc.disp(['k = ', num2str(min(p.k))], prms.verbalflag);
                dc.disp(['w = ', num2str(p.w)], prms.verbalflag);
                prms.plotNum = 1;
                figure(prms.plotNum); clf
                prms = dc.plot_timecourses(d, prms, p);
                if strcmp(lower(d.datatype), 'vep_kids')
                    text(1,.1, num2str(I.age));
                end
                disp(p)
                %  set(gcf, 'Position', [100 500 1300 300])
                set(gcf, 'Position', [750 500 1746 800])
                set(gcf, 'Position', [33         443        1369         241]);
                title([d.datatype, ' ', I.stimtype{1},' ', p.group{1}, ' ', p.id{1}, ' ', p.loc{1}], 'Interpreter','none');
                drawnow;
                if prms.savePlotOn
                    savefig(gcf,strcat(['data', filesep ,d.datatype, filesep, 'figures',  filesep, I.stimtype{1},'_',p.group{1}, '_',  p.id{1},'_', p.loc{1}]));
                    print(['data' filesep d.datatype filesep, 'figures', filesep, I.stimtype{1}, '_', p.group{1}, '_', p.id{1}, '_',p.loc{1}], '-dpng');
                    print(['data' filesep d.datatype filesep, 'figures', filesep, I.stimtype{1},'_', p.group{1}, '_', p.id{1}, '_',p.loc{1}], '-dmeta');
                end
            end
            % variance in the mean response, large numbers are better
            p.var = var([d.mean(1).response d.mean(2).response]);

            [err, errnoCost, d] = dc.fit_model( p, prms, d);
            p.corr = 1-errnoCost;
            dc.disp(p, prms.verbalflag);

            %% save all the data
            if ~exist('P', 'var')
                p.w = abs(p.w); % sometimes can be tiny neg numbers
                P = struct2table(p);
            else
                P(end+1, :) =  struct2table(p);
            end         
            if prms.pauseon; pause; end
        end % end of channels or ROIs]
    end % end of participants
    if exist('P')
        filename = ['data' filesep I.datatype{1},'_',prms.modelname, '_', I.stimtype{1}];
        if prms.averagehemis == true & size(Locs, 1)>1
            filename = strcat(filename, '_avehemi');
        end
        if prms.flipcontrast
            filename = strcat(filename, '_flipContrast');
        end
        writetable(P,strcat(filename, datetimestr, '.xlsx'));
    end
end
%cd(curDir);
%bold_correlations;
