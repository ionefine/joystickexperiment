
function plot_params_mri_ind()
close all
roi_list = {'V1', 'V2', 'V3'}; %,'vep_allchan',};
group_list = {'NS', 'BD' ,'AM'};
stim_list = {'congruent'};
param_list  ={'min_k', 'w'};

filename = ['data' filesep 'bold_',stim_list{1},'_model_fit_table.xlsx'];
A = readtable(filename);

for r = 1:3
    TMP = A(strcmp(A.loc, roi_list{r}), :);
    [r_kw(r), p_kw(r)] = corrcoef(TMP.min_k, TMP.w);
end

for p = 1:2
 V1 = A(strcmp(A.loc, roi_list{r}), :);


for g = 1:3
    for r = 1:3
        TT(g, r) = A(strcmp(A.group, group_list{g}) & strcmp(A.loc, roi_list{r}), :);
    end
end

color = colormap([1 0 0; .5 0 .75; 0 0 1]);

for r = 1:3

[r, p] = corrcoef(T.min_k, T.w)
p = polyfit(T.min_k,T.w,1);
plot(.7:.1:1, polyval(p, .7:.1:1)); hold on
for g = 1:length(group_list)
    G = T(strcmp(T.group, group_list{g}), :);
    id_List = unique(G.id)
    for i = 1:length(id_List)
        I = G(strcmp(G.id, id_List{i}), :);
        ph1 = scatter(I.min_k+0.01*(rand-.5), I.w+0.01*(rand-.5), 24, 'o'); hold on
        set(ph1,  'MarkerFaceColor',color(g,:), 'MarkerEdgeColor',color(g,:)); hold on
        ph1.MarkerFaceAlpha = .5;  ph1.MarkerEdgeAlpha = .5;
        xlabel('attenuation (k)');
        ylabel('binocular integration (w)');
        set(gca, 'XLim', [ 0 1]);
        set(gca, 'YLim', [ 0 1]);
    end
end



figure(2)
for p = 1:length(param_list)
    figure(p+1)
    for g = 1:length(group_list)
        subplot(1,3,1)
        ph1 = scatter(G1.min_k+0.01*(rand-.5), G2.min_k+0.01*(rand-.5), 24, 'o'); hold on
        set(ph1,  'MarkerFaceColor',color(g,:), 'MarkerEdgeColor',color(g,:)); hold on
        ph1.MarkerFaceAlpha = .5;  ph1.MarkerEdgeAlpha = .5;
        plot([0 1], [0 1], 'k')
        ph1.MarkerFaceAlpha = .5;  ph1.MarkerEdgeAlpha = .5;
        xlabel('V1'); ylabel('V2'); 
        set(gca, 'XLim', [-.1 1.1]); set(gca, 'YLim', [-.1 1.1]);axis square
         set(gca, 'XTick', 0:.25:1); set(gca, 'YTick', 0:.25:1);
        title('k')
        subplot(1,3,2)
        ph1 = scatter(G1.w+ 0.01*(rand-.5), G2.w+ 0.01*(rand-.5), 24, 'o'); hold on
        set(ph1,  'MarkerFaceColor',color(g,:), 'MarkerEdgeColor',color(g,:)); hold on
        plot([0 1], [0 1], 'k')
        ph1.MarkerFaceAlpha = .5;  ph1.MarkerEdgeAlpha = .5;
        xlabel('V1'); ylabel('V2'); 
        set(gca, 'XLim', [-.1 1.1]); set(gca, 'YLim', [-.1 1.1]); 
         set(gca, 'XTick', 0:.25:1); set(gca, 'YTick', 0:.25:1);
        axis square
        title('w'); set(gca, 'XTick', 0:.25:1);
    end
end

