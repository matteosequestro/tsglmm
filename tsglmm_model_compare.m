function modcompare_time_series_glmm(mods_to_compare, criterion)
% ===========================================================================
% Compare GLMM model fits across a time series
% 
% INPUTS:
%   ics       : matrix of model information criteria (nModels x nTimepoints)
%   criterion : string indicating the information criterion, e.g., 'BIC' or 'AIC'
%
% This function visualizes:
%   1. Average criterion difference from best model
%   2. Proportion of samples where each model is best
%   3. Which model is best across time
%   4. Delta between best and second-best models across time
% ===========================================================================

nModels = length(mods_to_compare);

% Convert criterion to char if it's a string
if isa(criterion, "string")
    criterion = convertStringsToChars(criterion);
end

% Extract required ICS
for md = 1 : nModels
    tIC = mods_to_compare(md).infcriteria.(criterion);
    if md == 1
        ics = NaN(nModels, length(tIC));
    end
    ics(md,:) = tIC;
end

% Unpack model names and fill if empty
modnames = [mods_to_compare.modname];
for ii = 1 : length(modnames)
    if strcmp(modnames(ii), "")
        modnames(ii) = "Mod" + ii;
    end
end

%% Compute average criterion differences
mean_ics = mean(ics, 2);              % mean across time for each model
delta_ics = min(mean_ics) - mean_ics; % difference from best model

% Identify best and second-best models
[~, idx_sorted] = sort(mean_ics); 
idx_two_best = idx_sorted(1:2);
delta_ics_winner = ics(idx_two_best(1),:) - ics(idx_two_best(2),:); % delta across time

%% Identify the best model at each timepoint
[~, bestmod] = min(ics, [], 1); % index of best model at each timepoint

% Compute proportion of times each model is best
prop_being_best = NaN(nModels, 1);
for md = 1 : nModels
    prop_being_best(md) = mean(bestmod == md);
end

%% Plotting
figure
tiledlayout(2,2); % 2x2 layout

% 1) Average criterion differences
nexttile(1);  % [rows columns] span
bar(delta_ics)
ylabel(['\Delta ', criterion, ' (from best)'])
xlabel('Model')
xticklabels(modnames)
title('Average criterion across time')
yline(0, 'LineStyle', '--')
ylim([min(delta_ics)-.1, .1])

% 2) Proportion of timepoints where each model is best
nexttile(2)
bar(prop_being_best)
ylabel({'Proportion of samples', 'being the best model'})
xlabel('Model')
xticklabels(modnames)
title('Proportion of samples being best')

% 3) Best model across time
nexttile(3)
hold on
for mm = 1 : nModels
    plot(bestmod == mm)
end
xlim([1, length(bestmod)])
ylim([-0.5, 1.5])
ylabel('Is this model the best?')
xlabel('Time (samples)')
legend(modnames, 'Location','northeast')
title('Best model across time')
hold off

% 4) Delta between best and second-best models across time
first_name  = convertStringsToChars(modnames(idx_two_best(1)));
second_name = convertStringsToChars(modnames(idx_two_best(2)));
nexttile(4)
plot(delta_ics_winner, 'LineWidth', 1.5)
ylabel(['\Delta ', criterion, ' (best - second best)'])
xlabel('Time (samples)')
yline(0, '--')
title(['\Delta', criterion, ' over time (', first_name, '-', second_name, ')'])
xlim([1, length(bestmod)])





end
