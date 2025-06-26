function f_dv_low_dim_proj_trial(app)

% app stuff
ops = app.ops;
params = f_dv_gather_params(app);

%tn_all = f_dv_get_trial_number(params);
[data, title_tag] = f_dv_get_data_by_mouse_selection(app.data, params);
%[region_num, reg_tag, leg_list] = f_dv_get_region_sel_val(params, ops);
plot_deets = 0;

if ~strcmpi(params.responsive_cells_select, 'All')
    params.responsive_cells_select = 'Resp marg';
end

num_dsets = size(data,1);

[tn0, plot_idx] = f_dv_get_trial_number(app, data.MMN_freq{1});
tn00 = tn0(1,:);
if strcmpi(data.paradigm{1}, 'FG_mmn')
    if sum(sum(tn00 == [1; 10])) == 2
        plot_idx = [plot_idx, {[find(tn00 == 1), find(tn00 == 10)]}];
    end
end

[num_gr, ~] = size(tn0);

data_all = cell(num_dsets,num_gr);
for n_dset = 1:num_dsets
    ddata = data(n_dset,:);
    stats1 = cat(1,ddata.stats{params.planes});
    tn1 = f_dv_get_trial_number(params, ddata.MMN_freq{1});
    
    for n_gr = 1:num_gr
        if 1%~sum(tn1(n_gr,:) == 0)
            [~, resp_vals] = f_dv_get_resp_vals_cells(stats1, tn1(n_gr,:), params);
            data_all{n_dset, n_gr} = cat(2,resp_vals{:});
        end
    end
end

data_all2 = cat(1,data_all{:});
data_all2 = f_dv_fix_nan_trials(data_all2, params.non_handle_method);

dred_params.method = params.dim_red_method;
dred_params.dist_metric = params.distance_method;
dred_params.subtract_mean = params.subtract_mean;
dred_params.scale_by_var = params.scale_by_var;
dred_params.plot_subtrat_mean = params.plot_subtracted_mean;
[lr_data, residual_var, residual_var_pca, subtr_mean] = f_dv_run_dred(data_all2, dred_params);

PR = (sum(residual_var_pca).^2)/(sum(residual_var_pca.^2));

title_tag1 = sprintf('%s; %s', title_tag, params.dim_red_method);
if strcmpi(params.method, 'isomap')
    title_tag1 = sprintf('%s; dist %s', title_tag1, params.distance_method);
end
title_tag2 = sprintf('%s; resp %s', title_tag1, params.responsive_cells_select);

if params.plot_subtracted_mean
    figure()
    plot(subtr_mean);
    xlabel('data axis');
    ylabel('mean magnitude');
    title(sprintf('subtracted mean; %s', title_tag2))
end

figure; hold on;
plot([0, 1:numel(residual_var_pca)], [1; residual_var_pca], 'o-', 'Linewidth', 2);
if ~strcmpi(params.method, 'pca')
    plot([0, 1:numel(residual_var)], [1; residual_var], 'o-', 'Linewidth', 2);
    legend('PCA', params.method)
end
title(sprintf('Residual variance proj trials; PR=%.2f; %s', PR, title_tag2), 'interpreter', 'none');
xlabel('number components used');
ylabel('Residual variance');

if strcmpi(params.num_axes_plot, '2')
    num_plots = ceil(params.num_comp_plot/2);
    % 2 comp  bs
    for n_pl = 1:num_plots
        pcs = [(n_pl-1)*2+1, (n_pl-1)*2+2];
        title_tag3 = sprintf('low rank proj trials; pl%d; %s', n_pl, title_tag2);

        f_dv_plot2_pc(lr_data, tn00, pcs, title_tag3, ops.context_types_all_colors2, plot_idx, params.render_painters);
        if strcmpi(params.method, 'pca')
            xlabel(sprintf('PC %d', pcs(1)));
            ylabel(sprintf('PC %d', pcs(2)));
        else
            xlabel(sprintf('Comp %d', pcs(1)));
            ylabel(sprintf('Comp %d', pcs(2)));
        end
    end
elseif strcmpi(params.num_axes_plot, '3')
    num_plots = ceil(params.num_comp_plot/3);
    % 3 comp  bs
    for n_pl = 1:num_plots
        pcs = [(n_pl-1)*3+1, (n_pl-1)*3+2, (n_pl-1)*3+3];
        title_tag3 = sprintf('low rank proj trials; pl%d; %s', n_pl, title_tag2);
        
        f_dv_plot3_pc(lr_data, tn00, pcs, title_tag3, ops.context_types_all_colors2, plot_idx, params.shadow_on3d, params.shadow_axis_locs, params.render_painters, params.grid_on, params.reverse_xyz)
    end
end

if plot_deets
    do_bar = 0;
    
    gcolor = gray(12);
    if do_bar
        y_max1 = round(max(lr_data(:)),1);
        y_min1 = round(min(lr_data(:)),1);
        num_bar = 3;
        figure()
        for n_pc = 1:num_bar
            subplot(num_bar,1,n_pc)
            bar(1:10, lr_data(:,n_pc))
            ylim([y_min1, y_max1])
        end
    else
        % plot trial variances within components
        figure(); hold on;
        for n_pc = 1:10
            plot(1:10, lr_data(:,n_pc), 'o-', color=gcolor(n_pc,:))
        end
        xlabel('trials')
        title('PCs magnitudes vs trials')
    end
    
    % plot trial variances within components
    figure(); hold on;
    for n_pc = 1:10
        plot(1:10, abs(lr_data(:,n_pc)), 'o-', color=gcolor(n_pc,:))
    end
    xlabel('trials')
    title('PCs abs magnitudes vs trials')

    figure(); hold on
    for n_tn = 1:10
        plot(1:10, (lr_data(n_tn,:)), 'o-', color=ops.context_types_all_colors2{tn_all(n_tn)})
    end
    title('Trial participation per PCs')
    xlabel('PCs')
    
    figure(); hold on
    for n_tn = 1:10
        plot(1:10, abs(lr_data(n_tn,:)), 'o-', color=ops.context_types_all_colors2{tn_all(n_tn)})
    end
    title('Trial abs participation per PCs')
    xlabel('PCs')
end
end