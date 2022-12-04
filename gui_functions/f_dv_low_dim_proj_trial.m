function f_dv_low_dim_proj_trial(app)

method = app.DimredmethodDropDown.Value;
dist_metric = app.LDdistmethodDropDown.Value; % pca isomap

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

num_dsets = size(data,1);

tn_all = f_dv_get_trial_number(app);

title_tag1 = sprintf('%s; %s', title_tag, method);
if strcmpi(method, 'isomap')
    title_tag1 = sprintf('%s; dist %s', title_tag1, dist_metric);
end

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

[lr_data2d, lr_data3d, residual_var, residual_var_pca] = f_dv_run_dred(data_all2, method, dist_metric);

title_tag2 = sprintf('%s; resp %s', title_tag1, resp_cell_sel);

figure; hold on;
plot([0, 1:numel(residual_var_pca)], [1; residual_var_pca], 'o-', 'Linewidth', 2);
if ~strcmpi(method, 'pca')
    plot([0, 1:numel(residual_var)], [1; residual_var], 'o-', 'Linewidth', 2);
    legend('PCA', method)
end
title(sprintf('Residual variance proj trials; %s', title_tag2), 'interpreter', 'none');
xlabel('number components used');
ylabel('Residual variance');

figure; hold on
plot(lr_data2d(:,1), lr_data2d(:,2), 'o-k', 'Linewidth', 1)
for n_tn = 1:numel(tn_all)
    plot(lr_data2d(n_tn, 1), lr_data2d(n_tn, 2), '.', 'color', app.ops.context_types_all_colors2{tn_all(n_tn)}, 'LineWidth', 2, 'MarkerSize', 20)
end
title(sprintf('low rank proj trials 2d; %s', title_tag2), 'interpreter', 'none');

title_tag3 = sprintf('low rank proj trials; %s', title_tag2);
f_dv_plot3_pc2(lr_data3d, tn_all, [], title_tag3, app.ops.context_types_all_colors2)

end