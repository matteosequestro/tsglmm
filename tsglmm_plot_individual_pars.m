function tsglmm_plot_individual_pars(modelout, parn, use_par_ci)


parnames_to_plot = modelout.pars.parnames;
parnames_to_plot = cellfun(@(x) replace(x, '_', ' '), parnames_to_plot, 'UniformOutput', false);

if nargin > 2; 
    use_par_ci = use_par_ci; 
else
    use_par_ci = 0; 
end

if nargin > 1

    indest = modelout.full_ind_estimates(:,:,parn);
    groupest = modelout.pars.estimates(parn,:);

    if use_par_ci
        lowerbound = modelout.pars.lower_ci(parn,:);
        upperbound = modelout.pars.upper_ci(parn,:);
    else
        tse = std(indest) ./ sqrt(height(indest));
        lowerbound = groupest - tse;
        upperbound = groupest + tse;
    end



    time = 1 : width(indest);

    figure
    plot(indest', 'Color', [0 0 0 .2])
    hold on
    patch([time, fliplr(time)], [lowerbound, fliplr(upperbound)], ...
        [.5,.5,.5], 'FaceAlpha',0.5, 'EdgeColor','none', 'HandleVisibility', 'off')
    plot(groupest, 'Color', [0 0 0], 'LineWidth', 2)
    yline(0, 'r--', 'LineWidth', 1.5)
else
    
    n_pars = size(modelout.full_ind_estimates, 3);

    nrows_tiles = ceil(sqrt(n_pars));
    ncols_tiles = ceil(n_pars / nrows_tiles);

    % Prepare the layout (keep it tight)
    figure
    tld             = tiledlayout(nrows_tiles, ncols_tiles);
    tld.TileSpacing = 'compact';
    tld.Padding     = 'compact';
    % Common labels
    ylabel(tld, 'parameter estimate (a.u.)');
    xlabel(tld, 'time (s)');

    for parn = 1 : n_pars
        indest = modelout.full_ind_estimates(:,:,parn);
        groupest = modelout.pars.estimates(parn,:);
        
        if use_par_ci
            lowerbound = modelout.pars.lower_ci(parn,:);
            upperbound = modelout.pars.upper_ci(parn,:);
        else
            tse = std(indest) ./ sqrt(height(indest));
            lowerbound = groupest - tse;
            upperbound = groupest + tse;
        end

        time = 1 : width(indest);
        nexttile
        plot(indest', 'Color', [0 0 0 .2])
        hold on
        patch([time, fliplr(time)], [lowerbound, fliplr(upperbound)], ...
            [.5,.5,.5], 'FaceAlpha',0.5, 'EdgeColor','none', 'HandleVisibility', 'off')
        plot(groupest, 'Color', [0 0 0], 'LineWidth', 2)
        yline(0, 'r--', 'LineWidth', 1.5)
        xlim([1, length(time)])

        title(parnames_to_plot{parn})
    end



end


