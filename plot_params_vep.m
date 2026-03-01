
datatype_list = {'vep'}; %,'vep_allchan',};
group_list = {'NS', 'BD' ,'AM'};
param_list  ={ 'err_noCost','min_k', 'w'}; % 
stim_list = {'orthogonal'};
% make sure err_noCost first
lim_list = [ 0 .13; -.5 1.1; -.5 1.1; log(.7) log(10); log(.7) log(10)];
color = colormap([1 0 0; .5 0 .75; 0 0 1]);

for dd = 1:length(datatype_list)
    for p = 1:length(param_list)
        figure(20+p); set(gcf, 'Name', [datatype_list{dd},' ', param_list{p}])
        for g = 1:length(group_list)
            for ss = 1:length(stim_list)
                filename = ['data' filesep datatype_list{dd},'_', stim_list{ss}, '_model_fit_table.xlsx'];
                T = readtable(filename);
                G = T(strcmp(T.group, group_list{g}), :);

                eval(['y = G.', param_list{p}, ';']);

                % remove bad fits
                if strcmp(param_list{p}, 'err_noCost') && length(y)>10
                    bvals = isoutlier(y);
                elseif ~exist('bvals')
                    bvals = zeros(size(y));
                end
                y(find(bvals)) = NaN;
                y(isoutlier(y)) = NaN; % remove outliers

                if contains(param_list(p), 'U')
                    y = log(y);
                    y(y<0) = 0.0001;
                end

                ph1 = scatter( g+.05*randn(size(y)), abs(y), 14,'o');
                set(ph1,  'MarkerFaceColor',color(g,:), 'MarkerEdgeColor',color(g,:)); hold on
                ph1.MarkerFaceAlpha = .5;  ph1.MarkerEdgeAlpha = .5;
                
                ph2 = scatter(g, nanmedian(abs(y)),50,'s');
                set(ph2, 'MarkerFaceColor',color(g,:), 'MarkerEdgeColor',color(g,:));
                ph2.MarkerFaceAlpha = .5;  ph2.MarkerEdgeAlpha = .3;
                line([g g], prctile(y, [25 75]),'Color', color(g,:), 'LineWidth',.5);
                
                set(gca, 'XLim', [.5 3.5]);
                set(gca, 'XTick', [1 2]);
                set(gca, 'XTickLabel', stim_list);
                set(gca, 'YLim', lim_list(p, :));

                if contains(param_list(p), 'U')
                    set(gca, 'YTick', log([0.1 1 10 ]));
                    dc.logy2raw()
                end
                % xlabel(stim_list{1})
                % ylabel(stim_list{2})
                % set(gca, 'XTick', 1:3, 'XTickLabel',group_list)
                title(group_list{g}, 'Interpreter','none')
                set(gcf, 'Position', [1892  849 560  420])
            end
        end
    end
end

