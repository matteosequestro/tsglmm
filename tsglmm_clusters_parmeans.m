function clust_par_means = tsglmm_clusters_parmeans(modelout)

% Compute mean parameter estimates for each individual and significant cluster
parnames = modelout.pars.parnames;
npar = length(parnames);

pars_series = modelout.full_ind_estimates;
alpha = modelout.alpha;

obs_clusters_sum = modelout.obs_clusters_sum;

if isempty(modelout.obs_clusters_sum)
    warning("no cluster dected")
else
    % Loop through parameters
    out_all_clust = []; % preallocate array
    variable_names = [];
    for par = 1 : npar

        % Unpack current parameter
        this_parameter      = pars_series(:,:,par);

        % Loop through significant clusters and compute mean pars
        this_par_clusters       = obs_clusters_sum(strcmp(obs_clusters_sum.pred, parnames{par}),:);
        significant_clusters    = find(this_par_clusters.prob < alpha);
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
            variable_names  = [variable_names, change_text([var_name, '_clst', num2str(tc)])];
        end
        % Compute cluster means for this parameter
        out_all_clust = [out_all_clust, out_par_means];
    end

    % Convert in nicer table for output
    clust_par_means    = array2table(out_all_clust, 'VariableNames', variable_names);
    clust_par_means.id = table2array(modelout.ids);
    clust_par_means    = clust_par_means(:, [end, 1:end-1]); % just because i like to have the id column as first
end
end