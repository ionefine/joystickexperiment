cd('C:\Users\Ione Fine\OneDrive - UW\Documents\code\amblyopia');
close all
roi_list = {'V1', 'V2', 'V3'}; %,'vep_allchan',};
group_list = {'NS', 'BD' ,'AM'};
stim_list = {'congruent'};
param_list  ={'err_noCost',  'w'};%

roiStr = {'V1', 'V2', 'V3'};
color = [.2  1 .2; 0 .7 .1; 0 .3 0];
allfilesP = dir(['*', 'psychophysics', '*','.xlsx']);
allfilesB = dir(['*', 'bold_weighted', '*','.xlsx']);

smallM = 30;
bigM = 60;
for ps = 1:length(allfilesP)
    for bd = 1:length(allfilesB)
        figure
        for r = 1:3
            P = readtable(allfilesP(ps).name);
            B = readtable(allfilesB(bd).name);
            P.w = abs(P.w);
            B.w = abs(B.w);
            
            subplot(1, 3, r)
            BR = B(strcmp(B.loc, roi_list{r}), :);
            % ind = find(BR.err_noCost<.035);
            % BR = BR(ind, :);
            % P = P(ind, :);

            % AMB
            amb = strcmp(BR.group, 'AM');
            x = P.w(amb);
            y = BR.w(amb);
            s_a=scatter(x, y, smallM+20, 'd'); hold on
            set(s_a, 'MarkerFaceColor', color(3,:));
            set(s_a, 'MarkerEdgeColor', color(3,:));
            set(s_a, 'MarkerFaceAlpha', .5);
            set(s_a, 'MarkerEdgeAlpha', .3);
                
            % BD
             bdr = strcmp(BR.group, 'BD');
            x = P.w(bdr);
            y = BR.w(bdr);
            s_bdr=scatter(x, y, smallM+20, 's'); hold on
            set(s_bdr, 'MarkerFaceColor', color(2,:));
            set(s_bdr, 'MarkerEdgeColor', color(2,:));
            set(s_bdr, 'MarkerFaceAlpha', .5);
            set(s_bdr, 'MarkerEdgeAlpha', .3);

            % NS
            ns = strcmp(BR.group, 'NS');
            x = P.w(ns);  mn_x_ns = mean(x);ste_x = std(x)./sqrt(length(x));
            y = BR.w(ns); mn_y_ns = mean(y);ste_y = std(y)./sqrt(length(y));
            xs = x+0.15*(rand(size(x))-.5)
            s_ns=scatter(xs, y, smallM);

            line([mn_x_ns mn_x_ns], [mn_y_ns-ste_y mn_y_ns+ste_y], 'Color',color(1,:));
            line([mn_x_ns-ste_x mn_x_ns+ste_x], [mn_y_ns mn_y_ns],'Color',color(1,:));

            set(s_ns, 'MarkerFaceColor', color(1,:));
            set(s_ns, 'MarkerEdgeColor', color(1,:));
            set(s_ns, 'MarkerFaceAlpha', .5);
            set(s_ns, 'MarkerEdgeAlpha', .3);
            
            [rho,pval] = corr(P.w(ns), BR.w(ns));
            disp(['Corr roi = ', roi_list{r}, ' p = ', num2str(pval)]);
            [rho,pval] = corr(P.w(ns), BR.w(ns), 'Type', "Spearman");
            disp(['Spearman roi = ', roi_list{r}, ' p = ', num2str(pval)]);
            xlabel('psychophysics')
            ylabel('BOLD')
            plot([0 1], [0 1], 'k--'); hold on

            s_ns_mn=scatter(mn_x_ns, mn_y_ns, bigM); hold on
            set(s_ns_mn, 'MarkerFaceColor', color(1,:));
            set(s_ns_mn, 'MarkerEdgeColor', color(1,:));
            set(s_ns_mn, 'MarkerFaceAlpha', 1);
            set(s_ns_mn, 'MarkerEdgeAlpha', 1);
            set(gca, 'XLim', [-.1 1.1]);
            set(gca, 'YLim', [-.1 1.1]);
            axis square
            set(gcf, 'Position', [680   288   951   710])

        end
    end
end
cd('C:\Users\Ione Fine\OneDrive - UW\Documents\code\amblyopia');