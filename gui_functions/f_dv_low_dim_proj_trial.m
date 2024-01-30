function f_dv_low_dim_proj_trial(app)

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

num_dsets = size(data,1);

tn_all = f_dv_get_trial_number(app);



data_all = cell(num_dsets,1);
for n_dset = 1:num_dsets
    stats1 = cat(1,data(n_dset,:).stats{n_pl});

    if strcmpi(app.ResponsivecellsselectDropDown.Value, 'All')
        resp_cell_sel = 'All';
    else
        resp_cell_sel = 'Resp marg';
    end

    [~, resp_vals] = f_dv_get_resp_vals_cells(app, stats1, tn_all, [], resp_cell_sel);

    data_all{n_dset} = cat(2,resp_vals{:});
end

data_all2 = cat(1,data_all{:});
hasnan1 = logical(sum(isnan(data_all2),2));
data_all2 = data_all2(~hasnan1,:);

params.method = app.DimredmethodDropDown.Value;
params.dist_metric = app.DistmethodDropDown.Value;
params.subtract_mean = app.subtractmeanCheckBox.Value;
params.scale_by_var = app.scalebyvarCheckBox.Value;
params.plot_subtrat_mean = app.plotsubmeanCheckBox.Value;
[lr_data, residual_var, residual_var_pca, subtr_mean] = f_dv_run_dred(data_all2, params);

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
title(sprintf('Residual variance proj trials; %s', title_tag2), 'interpreter', 'none');
xlabel('number components used');
ylabel('Residual variance');

if strcmpi(app.NumplotaxesDropDown.Value, '2')
    num_plots = ceil(app.numcompplotSpinner.Value/2);
    lim_max = 0;
    % 2 comp  bs
    for n_pl = 1:num_plots
        ax1 = (n_pl-1)*2+1;
        ax2 = (n_pl-1)*2+2;
        lim_max = max([ceil(max([max(lr_data(:,ax1)) - min(lr_data(:,ax1)), max(lr_data(:,ax2)) - min(lr_data(:,ax2))])*1.05*100)/100, lim_max]);
        cent1 = (max(lr_data(:,ax1)) + min(lr_data(:,ax1)))/2;
        cent2 = (max(lr_data(:,ax2)) + min(lr_data(:,ax2)))/2;
        figure; hold on
        if app.plotaxesCheckBox.Value
            plot(linspace(cent1-lim_max/2, cent1+lim_max/2,10), zeros(10,1), color=[0.7, 0.7, 0.7], LineWidth=1);
            plot(zeros(10,1), linspace(cent2-lim_max/2, cent2+lim_max/2,10), color=[0.7, 0.7, 0.7], LineWidth=1);
        end
        pl1 = plot(lr_data(:,ax1), lr_data(:,ax2), 'o-k', 'Linewidth', 1);
        for n_tn = 1:numel(tn_all)
            plot(lr_data(n_tn, ax1), lr_data(n_tn, ax2), '.', 'color', app.ops.context_types_all_colors2{tn_all(n_tn)}, 'LineWidth', 2, 'MarkerSize', 20)
        end
        xlabel(sprintf('PC%d', ax1));
        ylabel(sprintf('PC%d', ax2));
        title(sprintf('low rank proj trials 2d; pl%d; %s', n_pl, title_tag2), 'interpreter', 'none');
        if app.equalizeaxesCheckBox.Value
            pl1.Parent.XLim = [cent1-lim_max/2, cent1+lim_max/2];
            pl1.Parent.YLim = [cent2-lim_max/2, cent2+lim_max/2];
        end
    end
elseif strcmpi(app.NumplotaxesDropDown.Value, '3')
    num_plots = ceil(app.numcompplotSpinner.Value/3);
    % 3 comp  bs
    for n_pl = 1:num_plots
        pcs = [(n_pl-1)*3+1, (n_pl-1)*3+2, (n_pl-1)*3+3];
        title_tag3 = sprintf('low rank proj trials; pl%d; %s', n_pl, title_tag2);
        f_dv_plot3_pc2(lr_data, tn_all, pcs, title_tag3, app.ops.context_types_all_colors2)
    end
end

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