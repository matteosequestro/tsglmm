




function [out, f1] = tsgmm_plot_clust_par_means(data)

parnames = data.Properties.VariableNames;



data = table2array(data ); % example data
means = mean(data);
ses = std(data) ./ sqrt(height(data));



f1 = figure; 
hold on
for i = 1:width(data)
    x = i + 0.1*(rand(size(data,1),1)-0.5); % jitter
    scatter(x, data(:,i), 15, ...
        'k', 'filled', ...
        'MarkerFaceAlpha', 0.5);
end

errorbar(1:width(data), means, ses, 'o',...
    'MarkerSize', 3, 'MarkerFaceColor', 'r', 'Color', 'r', ...
    'LineWidth', 1.2, 'CapSize', 0);
xlim([0.5 width(data)+.5])
set(gca, 'XTick', 1:width(data))
set(gca, 'XTicklabels', parnames )
box on
yline(0, '--')
ylabel('parameter estimate (a.u.)')



BF10 = NaN(width(data), 1);
for ii = 1 : width(data)
    BF10(ii) = bf.ttest(data(:,ii),0);
end

out.parnames = parnames;
out.BF10 = BF10;
out.BF01 = 1/BF10;
out.means = means;
out.ses = ses;




end