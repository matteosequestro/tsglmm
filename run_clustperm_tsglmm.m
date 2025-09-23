

function [modelout, obs_clusters_sum, outputText, outputStats, clust_par_means] = run_clustperm_glmm(data, formula, cfg)
% ----------------------------------------------------------------------------------------------------------------------------------------------
% Performs cluster based correction on parameter estimates from linear mixed
% effect model on a time-series
%
% INPUTS:
% data                  : full dataset having at least a column for ID,
%                           one column for each predictor in your model,
%                           and one column with the time series in a cell.
%                           Each row is a trial.
% formula               : the model formula in Wilkinson notation (i.e., R-like
%                           notation, e.g., 'y ~ 1 + x1 * x2 + (1 + x1 * x2 | id)'.
%                           It must be a string.
% cfg.niter             : number of permutations.
% cfg.wantplot_perm     : 1 if you want to plot time-series for each
%                           parameter with significant (blue) and
%                           non-significant (gray) clusters as shaded
%                           reptangles. 0 if you don't want such a plot.
% cfg.perm_alpha        : permutation alpha. Montecarlo p-values below such
%                           a value will be considered as significant.
% cfg.want_parallel_fit : 1 will fi the models across different time-points
%                           in parallel (much faster but computationally
%                           more expensive); 0 will do one point at the
%                           time.
% OUTPUTS:
% modelout              : ouput of the model on true dataset (i.e., not
%                           permuted. You will get a "pars" field
%                           (parameter estimate, lower and upper confidence
%                           interval, tvalues), a "criteria" field with
%                           log-likelihood of the model, as well as AIC,
%                           BIC and deviance criteria (for each sample
%                           point; and an "rsqrd" field with ordinary and
%                           adjusted R-squared
% obs_clusters_sum      : output of the cluster based permutation. For each
%                           coefficient in the model and each cluster
%                           detected on the true dataset (i.e., not
%                           permutedyou), get a monte-carlo
%                           p-value, first sample in the cluster, length
%                           of the cluster and t-mass.
% ----------------------------------------------------------------------------------------------------------------------------------------------


% Set defaults if not set in the cfg
plot_non_significant_clusters   = getOrDefault(cfg, 'plot_non_significant_clusters', 0);
glm_likelihood                  = getOrDefault(cfg, 'glm_likelihood', 'Gaussian');
cut_short_clusters              = getOrDefault(cfg, 'cut_short_clusters', 0);
minlength                       = getOrDefault(cfg, 'minlength', 0);
wantplot_perm                   = getOrDefault(cfg, 'wantplot_perm', 1);
modname                         = getOrDefault(cfg, 'modname', ""); % just for plotting purposes, e.g., if you need model comparison

% Take time
tic

% Unpack the fieldtrip cfg with something uninteresting
fieldtrip_cfg            = cfg.fieldtrip_cfg;
fieldtrip_cfg.uvar       = 1;  % subject
fieldtrip_cfg.ivar       = 2;  % condition
fieldtrip_cfg.spmversion = 'spm12';
fieldtrip_cfg.feedback   = 'no';

% Get defaults for fieldtrip
fieldtrip_cfg.numrandomization  = getOrDefault(fieldtrip_cfg, 'numrandomization', 1e3);              % number of permutations
fieldtrip_cfg.method            = getOrDefault(fieldtrip_cfg, 'method',           'montecarlo');     % use monte carlo sampling for 'p-values'
fieldtrip_cfg.statistic         = getOrDefault(fieldtrip_cfg, 'statistic',        'depsamplesT');    % paired t-test (equivalent to t-test agains 0 if contrasted with a null matrix of zeros)
fieldtrip_cfg.correctm          = getOrDefault(fieldtrip_cfg, 'correctm',         'cluster');        % correction method
fieldtrip_cfg.clusterstatistic  = getOrDefault(fieldtrip_cfg, 'clusterstatistic', 'maxsum');         % statistic to compute on permutations
fieldtrip_cfg.tail              = getOrDefault(fieldtrip_cfg, 'tail',             0);                % two-sided test to find clusters
fieldtrip_cfg.clustertail       = getOrDefault(fieldtrip_cfg, 'clustertail',      0);                % compare observed clusters against positive and negative permutation clusters
fieldtrip_cfg.clusteralpha      = getOrDefault(fieldtrip_cfg, 'clusteralpha',     0.05);             % alpha for cluster detection
fieldtrip_cfg.alpha             = getOrDefault(fieldtrip_cfg, 'alpha',            0.05);             % alpha for monte carlo testing
fieldtrip_cfg.neighbours        = getOrDefault(fieldtrip_cfg, 'neighbours',       []);               % no spatial neighbours since it's unidimensional

% Fit the model and find clusters in the observed data
disp('_________________________________________________')
disp('fitting the model...')
[modelout]          = time_series_glmm(data, formula, cfg);
fprintf('fitting done in %.2fs!\n', round(toc,2));

npar                = height(modelout.pars.estimates); % number of parameters (a.k.a., coefficients)
tslen               = width(modelout.pars.estimates); % length of the time-series
parnames            = modelout.pars.parnames;
modelout.formula    = formula;
modelout.modname    = modname;

% Unpack individual parameter estimates
pars_series = modelout.full_ind_estimates;

% Run Fieldtrip cluster based permutationn
n_subjs = height(pars_series);
timevec = 1 : width(pars_series);

% Design matrix for fieldtrip
design(1,:) = repmat(1 : n_subjs, 1, 2);
design(2,:) = [ones(1, n_subjs), ones(1, n_subjs) .* 2];
fieldtrip_cfg.design = design;

% Preallocate output arrays
outputText = cell(npar, 1);
outputStats = cell(npar, 1);

% Loop through parameters
disp('permuting...')
for pp = 1 : size(pars_series, 3)
    mySignal = pars_series(:,:,pp);

    % Shape parameters in fieldtrip format
    subjdata = cell(1, n_subjs);
    for subj = 1 : n_subjs
        subjdata{subj}.avg    = mySignal(subj,:);   % 1 x nTime of signal
        subjdata{subj}.time   = timevec;            % 1 x nTime 1 to last moment
        subjdata{subj}.label  = {'signal'};         % single channel label
        subjdata{subj}.dimord = 'chan_time';
    end

    % Create contrast data (only once, then it's the same for all pars)
    if pp == 1
        zero_data = subjdata;
        for dd = 1 : length(zero_data)
            zero_data{dd}.avg = zero_data{dd}.avg .* 0;
        end
    end

    % Run cluster based permutation (Maris and Oostenveld, 2007, Journal of Neuroscience Methods)
    outputText{pp} = evalc('stat = ft_timelockstatistics(fieldtrip_cfg, subjdata{:}, zero_data{:});');
    outputStats{pp} = stat;

    % Extract positive cluster info
    n_pclusters = length(stat.posclusters);
    for cl = 1 : n_pclusters
        clus_pos = find(stat.posclusterslabelmat == cl);
        stat.posclusters(cl).first      = clus_pos(1);
        stat.posclusters(cl).last       = clus_pos(end);
        stat.posclusters(cl).length     = length(clus_pos);
    end
    
    % Extract negative cluster info
    n_nclusters = length(stat.negclusters);
    for cl = 1 : n_nclusters
        clus_neg                        = find(stat.negclusterslabelmat == cl);
        stat.negclusters(cl).first      = clus_neg(1);
        stat.negclusters(cl).last       = clus_neg(end);
        stat.negclusters(cl).length     = length(clus_neg);
    end
    
    % Convert clusters to tables once
    posT = struct2table(stat.posclusters);
    negT = struct2table(stat.negclusters);

    % Merge positive and negative clusters
    if ~isempty(posT) && ~isempty(negT)
        clusters = [posT; negT];
    elseif ~isempty(posT)
        clusters = posT;
    elseif ~isempty(negT)
        clusters = negT;
    else
        continue  % skips to next iteration of the for-loop
    end
    clear posT negT

    clusters            = sortrows(clusters, 'first');                      % Sort clusters by timing
    clusters.clustern   = (1 : height(clusters))';                          % Add cluster number
    clusters.pred       = repmat(switchText(parnames{pp}), height(clusters), 1);        % Add predictor name

    % Export to the main dataset (or create it)
    if pp == 1
        obs_clusters_sum = clusters;
    else
        obs_clusters_sum = [obs_clusters_sum; clusters];
    end
end

% Cut short clusters
if cut_short_clusters
    obs_clusters_sum(obs_clusters_sum.length < minlength,:) = [];
end

modelout.obs_clusters_sum = obs_clusters_sum;

%% Plots
if wantplot_perm
  

    
end % cfg.wantplot_perm


%% Compute mean parameter estimates for each individual and significant cluster
if nargout > 4 % Do this only if it's required by the output to save time
    if isempty(obs_clusters_sum)
        clust_par_means = [];
    else
        % Loop through parameters
        out_all_clust = []; % preallocate array
        variable_names = [];
        for par = 1 : npar

            % Unpack current parameter
            this_parameter      = pars_series(:,:,par);

            % Loop through significant clusters and compute mean pars
            this_par_clusters       = obs_clusters_sum(strcmp(obs_clusters_sum.pred, parnames{par}),:);
            significant_clusters    = find(this_par_clusters.prob < fieldtrip_cfg.alpha);
            out_par_means           = NaN(height(this_parameter), length(significant_clusters));
            for tc = 1 : length(significant_clusters)
                this_cluster            = significant_clusters(tc);
                clFirst                 = this_par_clusters.first(this_cluster);
                clLength                = this_par_clusters.length(this_cluster);
                end_cluster = clFirst  + clLength - 1;
                if end_cluster > width(this_parameter)
                    end_cluster = width(this_parameter);
                end
                out_par_means(:,tc)     = mean(this_parameter(:, clFirst : end_cluster), 2);

                % Make parameters name suitable to be put as column names in
                % the output and add the cluster number
                var_name        = erase(parnames{par}, ['(', ')']);
                var_name        = replace(var_name, ':', 'X');
                variable_names  = [variable_names, switchText([var_name, '_clst', num2str(significant_clusters(tc))])];
            end
            % Compute cluster means for this parameter
            out_all_clust = [out_all_clust, out_par_means];
        end

        % Convert in nicer table for output
        clust_par_means    = array2table(out_all_clust, 'VariableNames', variable_names);
        clust_par_means.id = modelout.ids;
        clust_par_means    = clust_par_means(:, [end, 1:end-1]); % just because i like to have the id column as first
    end
end

% Display total elapsed time
fprintf('done in %.2fs!\n', round(toc,2));
disp('_________________________________________________')


end
