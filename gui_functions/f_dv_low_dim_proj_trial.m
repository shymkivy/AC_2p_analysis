function f_dv_low_dim_proj_trial(app)

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

plot_deets = 0;

shadow_axis_locs = [app.FlipshadowXCheckBox.Value, app.FlipshadowYCheckBox.Value, app.FlipshadowZCheckBox.Value] + 1;
reverse_xyz = [app.ReverseXCheckBox.Value, app.ReverseYCheckBox.Value, app.ReverseZCheckBox.Value];

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

if strcmpi(app.ResponsivecellsselectDropDown.Value, 'All')
    resp_cell_sel = 'All';
else
    resp_cell_sel = 'Resp marg';
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
    stats1 = cat(1,ddata.stats{n_pl});
    tn1 = f_dv_get_trial_number(app, ddata.MMN_freq{1});
    
    for n_gr = 1:num_gr
        if 1%~sum(tn1(n_gr,:) == 0)
            [~, resp_vals] = f_dv_get_resp_vals_cells(app, stats1, tn1(n_gr,:), [], resp_cell_sel);
            data_all{n_dset, n_gr} = cat(2,resp_vals{:});
        end
    end
end

data_all2 = cat(1,data_all{:});
data_all2 = f_dv_fix_nan_trials(data_all2, app.nanhandlemetDropDown.Value);

params.method = app.DimredmethodDropDown.Value;
params.dist_metric = app.DistmethodDropDown.Value;
params.subtract_mean = app.subtractmeanCheckBox.Value;
params.scale_by_var = app.scalebyvarCheckBox.Value;
params.plot_subtrat_mean = app.plotsubmeanCheckBox.Value;
[lr_data, residual_var, residual_var_pca, subtr_mean] = f_dv_run_dred(data_all2, params);

PR = (sum(residual_var_pca).^2)/(sum(residual_var_pca.^2));

title_tag1 = sprintf('%s; %s', title_tag, params.method);
if strcmpi(params.method, 'isomap')
    title_tag1 = sprintf('%s; dist %s', title_tag1, params.dist_metric);
end
title_tag2 = sprintf('%s; resp %s', title_tag1, resp_cell_sel);

if app.plotsubmeanCheckBox.Value
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

if strcmpi(app.NumplotaxesDropDown.Value, '2')
    num_plots = ceil(app.numcompplotSpinner.Value/2);
    % 2 comp  bs
    for n_pl = 1:num_plots
        pcs = [(n_pl-1)*2+1, (n_pl-1)*2+2];
        title_tag3 = sprintf('low rank proj trials; pl%d; %s', n_pl, title_tag2);

        f_dv_plot2_pc(lr_data, tn00, pcs, title_tag3, app.ops.context_types_all_colors2, plot_idx, app.FigrenderpaintersCheckBox.Value)
    end
elseif strcmpi(app.NumplotaxesDropDown.Value, '3')
    num_plots = ceil(app.numcompplotSpinner.Value/3);
    % 3 comp  bs
    for n_pl = 1:num_plots
        pcs = [(n_pl-1)*3+1, (n_pl-1)*3+2, (n_pl-1)*3+3];
        title_tag3 = sprintf('low rank proj trials; pl%d; %s', n_pl, title_tag2);
        
        f_dv_plot3_pc(lr_data, tn00, pcs, title_tag3, app.ops.context_types_all_colors2, plot_idx, app.shadowon3dCheckBox.Value, shadow_axis_locs, app.FigrenderpaintersCheckBox.Value, app.gridon3dCheckBox.Value, reverse_xyz)
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
        plot(1:10, (lr_data(n_tn,:)), 'o-', color=app.ops.context_types_all_colors2{tn_all(n_tn)})
    end
    title('Trial participation per PCs')
    xlabel('PCs')
    
    figure(); hold on
    for n_tn = 1:10
        plot(1:10, abs(lr_data(n_tn,:)), 'o-', color=app.ops.context_types_all_colors2{tn_all(n_tn)})
    end
    title('Trial abs participation per PCs')
    xlabel('PCs')
end
end