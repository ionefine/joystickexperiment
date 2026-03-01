function plot_params_mri()

roi_list = {'V1', 'V2', 'V3'}; %,'vep_allchan',};
group_list = {'NS', 'BD' ,'AM'};
stim_list = {'congruent'};
param_list  ={'err_noCost',  'w'};% 
% make sure err_noCost first
lim_list = [ 0 .08; -.5 1.1; -.5 1.1; log10(.7) log10(10); log10(.7) log10(10)];
color = [.2  1 .2; 0 .7 .1; 0 .3 0];

filename = ['data' filesep 'bold_weighted_binocular_norm_congruent_avehemi2601292349.xlsx']
T = readtable(filename);
loc_List = unique(T.loc);
for g = 1:length(group_list)
    clear y all_y gvals
    G = T(strcmp(T.group, group_list{g}), :);

    for p = 1:length(param_list)
        figure(p+10); set(gcf, 'Name', ['MRI ', param_list{p}])
        subplot(1,3, g);
        for c = 1:length(loc_List) % rois
            C = G(strcmp(G.loc, loc_List{c}), :);
            eval(['y = C.', param_list{p}, ';']);
            % remove bad fits
            if strcmp(param_list{p}, 'err_noCost') && length(y)>10
                gvals = ~isoutlier(y);
            elseif ~exist('gvals')
                gvals = ones(size(y));
            end
            y = y(find(gvals));

            y= y(find(~isoutlier(y))); % remove outliers

            if contains(param_list(p), 'U')
                y = log10(y);
                y(y<0) = 0.0001;
            end

            ph1 = scatter(c+.05*randn(size(y)), abs(y), 14,'o'); hold on
            set(ph1,  'MarkerFaceColor',color(g,:), 'MarkerEdgeColor',color(g,:)); hold on
            ph1.MarkerFaceAlpha = .5;  ph1.MarkerEdgeAlpha = .5;
            if g==1
                ph2 = scatter(c, median(abs(y)),50,'s');
            set(ph2, 'MarkerFaceColor',color(g,:), 'MarkerEdgeColor',color(g,:));
            ph2.MarkerFaceAlpha = .5;  ph2.MarkerEdgeAlpha = .3;
            line([c c], prctile(y, [25 75]),'Color', color(g,:), 'LineWidth',.5);
            end
            set(gca, 'YLim', lim_list(p, :))
            if contains(param_list(p), 'U')
                set(gca, 'YTick', log10([0.001 .1 1 10 100]));
                dc.logy2raw(10)
            end

            set(gca, 'XLim', [.5 3.5] );
            set(gca, 'XTick', 1:3, 'XTickLabel',roi_list);
            title(group_list{g})
            set(gcf, 'Position', [46   689   785   420])
        end
    end
end
% end
% figure(10)
% errorbar([1 2 3], [1 1.328627038 1.463394159], [0 0.070305466 0.110516075])
