

function output = tsglmm_diagnostics(modelout)


%%% Compute everything you need before plotting

% number of bootstraps to compute the confidence interval for posterior
% predictive check (default at 100, which is not ideal in general for
% bootstrapped CIs, but likely you don't even
if nargin < 2; nboots = 100; end

% Unpack VIF table
viftable    = modelout.vifs;
Predictor   = categorical(cellfun(@(x) replace(x, ':', 'X'), viftable{:,1}, 'UniformOutput', false));
vif         = viftable{:,2};         

% Unpack real and predicted values, compute correlations and CIs
y       = modelout.predictions.y;                                                                   % observed by-subject averages
yhat    = modelout.predictions.yhat;                                                                % predicted by-subject averages
r       = arrayfun(@(i) corr(y(:,i), yhat(:,i)), 1:size(y,2));                                      % by-sample correlation between true and predicted averages
ci      = arrayfun(@(i) corr_boot(y(:,i), yhat(:,i), nboots), 1:size(y,2), 'UniformOutput',false);  % bootstrapped correlations for the above mentioned
lower   = cellfun(@(x) x(1), ci);                                                                   % lower bound of confidence interval
upper   = cellfun(@(x) x(2), ci);                                                                   % upper confidence interval

% Unpack residuals
residuals = modelout.residuals;
stdz_residuals = modelout.stdz_residuals;
pooled = stdz_residuals(:);

% Extract some useful variables
x = 1:width(y);                                                                               % time vector to put as x axis

% Create the main layout (each tile is a criterion)
figure
tld = tiledlayout(3,2);
tld.TileSpacing = 'compact';
tld.Padding     = 'compact';
title(tld, replace(modelout.formula, '_', ' ')); % put the main formula as title

%% R-squared
nexttile
plot(modelout.rsqrd.ordinary)
hold on
plot(modelout.rsqrd.adjusted)
legend(["ordinary", "adjusted"])
xlabel("time (sample)")
ylabel("R²")
title("Goodness of fit")
xlim([1, max(x)])
ylim([min(min(modelout.rsqrd.ordinary), min(modelout.rsqrd.adjusted)), ...
    max(max(modelout.rsqrd.ordinary), max(modelout.rsqrd.adjusted))*1.1])
box on
hold off


%% Posterior predictive check (PPC)
nexttile
meany = mean(y);
meanyhat = mean(yhat);
plot(meany)
hold on
plot(meanyhat)
legend(["real grand average", "predicted grand average"], ...
    'Location', 'southeast')
ylabel("y")
xlabel("Time (sample)")
title("Posterior predictive check 1")
xlim([1, max(x)])
box on
ylim([min(min(meany), min(meanyhat)), ...
    max(max(meany), max(meanyhat))*1.1])

%%% Plot correlation + bootstrapped CI for predicted vs real across time
nexttile
% Shaded CI area
fill([x fliplr(x)], [lower fliplr(upper)], ...
    [0.8 0.8 0.8], 'EdgeColor','none', 'FaceAlpha',0.3);
hold on;

% Plot correlations
plot(x, r, 'k', 'LineWidth', 1.5);
legend(["95% bootstrapped confidence interval", "pearson r(real,predicted)"], ...
    'Location', 'southeast')
xlabel('Time (sample)');
ylabel('Correlation (r)');
ylim([min(lower-0.05), max(upper+0.05)]);
xlim([1, max(x)])
box on
title("Posterior predictive check 2")

%% Variance Inflation factors
% Draw rectangles first (background shading)
nexttile
hold on
barWidth = length(Predictor) ;
scatter(Predictor, vif, 40, 'k', 'filled')
if max(vif) >= 5
    rectangle('Position', [0.5, 0, barWidth, 5], ...
        "FaceColor", [0.5 0.7 0], 'FaceAlpha', 0.3, 'EdgeColor', "none")
    if max(vif)>= 10
        rectangle('Position', [0.5, 5, barWidth, 5], ...
            "FaceColor", [0.9 0.9 0], 'FaceAlpha', 0.3, 'EdgeColor', "none")
        rectangle('Position', [0.5, 10, barWidth,max(vif)-9], ...
            "FaceColor",[0.9 0 0], 'FaceAlpha', 0.3, 'EdgeColor', "none")
    else
        rectangle('Position', [0.5, 5, barWidth,  max(vif)-4], ...
            "FaceColor", [0.9 0.9 0], 'FaceAlpha', 0.3, 'EdgeColor', "none")
    end
else
    rectangle('Position', [0.5, 0, barWidth, max(vif)+1], ...
        "FaceColor", [0.5 0.7 0], 'FaceAlpha', 0.3, 'EdgeColor', "none")
end
ylabel('VIF')
xlabel('Predictor')
title("Variance Inflation Factors")
box on
hold off

%% Residuals
nexttile
histogram(pooled, 'NumBins', 50, 'Normalization', 'pdf')  % normalize
hold on
[f, xi] = ksdensity(pooled);
plot(xi, f, 'r--', 'LineWidth', 1)
xlabel('Standardized residual')
ylabel('Probability density')
title('Histogram of residuals')

nexttile
qqplot(pooled)
title('QQ-plot of residuals')


%% Create table for export
output.ppc_corrs        = table(r',lower',upper', 'VariableNames', ["r", "lowerci", "upperci"]);
output.rsqrd            = modelout.rsqrd;
output.vif              = viftable;
output.residuals        = modelout.residuals;
output.stdz_residuals   = modelout.stdz_residuals;

end



