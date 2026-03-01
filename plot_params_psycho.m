clear

datatype_list = {'psychophysics'}; %{'psychophysics'}; %,'vep_allchan',};
group_list = {'NS', 'BD' ,'AM'};
stim_list = {'congruent'};
param_list  ={'err_noCost', 'min_k', 'w'}; %
% make sure err_noCost first
lim_list = [ 0 .13; -.5 1.1; -.5 1.1; log10(.7) log10(10); log10(.7) log10(10)];
color = [.2  1 .2; 0 .7 .1; 0 .3 0];

for dd = 1:length(datatype_list)
   %filename = ['data' filesep 'bold_psychophysics_weighted_binocular_norm_congruent.xlsx '];
    filename = 'bold_psychophysics_weighted_binocular_norm_congruent2601292246.xlsx'
    T = readtable(filename);
    chan_list = unique(T.loc);
    C = T(strcmp(T.loc, chan_list{1}), :);
    for g = 1:length(group_list)
        G = C(strcmp(C.group, group_list{g}), :);
        clear gvals
        for p = 1:length(param_list)
            figure(p); set(gcf, 'Name', ['PSYCHO ', param_list{p}]);
            eval(['y = G.', param_list{p}, ';']);
            y = abs(y);
            % remove bad fits
            if strcmp(param_list{p}, 'err_noCost') && length(y)>10
                gvals = ~isoutlier(y);
            elseif ~exist('gvals')
                gvals = ones(size(y));
            end
        %    y = y(find(gvals));

       %     y= y(find(~isoutlier(y))); % remove outliers

            if contains(param_list(p), 'U')
                y = log10(y);
                y(y<0) = 0.0001;
            end

            ph1 = scatter(g+.03*randn(size(y)), y, 14,'o');
            set(ph1,  'MarkerFaceColor',color(g,:), 'MarkerEdgeColor',color(g,:)); hold on
            ph1.MarkerFaceAlpha = .5;  ph1.MarkerEdgeAlpha = .5;
          if g==1
              ph2 = scatter(g, median(y),50,'s');
            set(ph2, 'MarkerFaceColor',color(g,:), 'MarkerEdgeColor',color(g,:));
            ph2.MarkerFaceAlpha = .5;  ph2.MarkerEdgeAlpha = .3;
            ph4 = scatter(g, prctile(y, [50]), 18, 's');
            if p==3 && dd==1
                prctile(y, [50])
            end
            set(ph4, 'MarkerFaceColor', color(g,:), 'MarkerEdgeColor',color(g,:)); hold on
            line([g g], prctile(y, [25 75]),'Color', color(g,:), 'LineWidth',.5);
          end 
          set(gca, 'YLim', lim_list(p, :));
            if contains(param_list(p), 'U')
                set(gca, 'YTick', log10([0.1 1 10 ]));
                dc.logy2raw(10)
            end
            set(gca, 'XLim', [.5 3.25] );
            set(gca, 'XTick', 1:3, 'XTickLabel',group_list)
            title(datatype_list{dd}, 'Interpreter','none')
             set(gcf, 'Position', [46   689   785   420])

        end
    end
end
