function f_dv_low_dim_proj_cell(app)

method = app.DimredmethodDropDown.Value;
dist_metric = app.DistmethodDropDown.Value; % pca isomap

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

cell_labels = cell(num_dsets,1);
data_all = cell(num_dsets,1);
for n_dset = 1:num_dsets
    stats1 = cat(1,data(n_dset,:).stats{n_pl});

    if strcmpi(app.ResponsivecellsselectDropDown.Value, 'All')
        resp_cell_sel = 'All';
    else
        resp_cell_sel = 'Resp marg';
    end

    [sel_cells, resp_vals, resp_vals2, ~, resp_cells] = f_dv_get_resp_vals_cells(app, stats1, tn_all, [], resp_cell_sel);
    
    sel_idx = logical(sum(sel_cells,2));
    resp_idx = logical(sum(resp_cells,2));
    [~, max_idx] = max(resp_vals2, [], 2);
    
    max_idx(~resp_idx) = 0;
    
    if strcmpi(app.ResponsivecellsselectDropDown.Value, 'All')
        cell_labels{n_dset} = max_idx;
    else
        cell_labels{n_dset} = max_idx(sel_idx);
    end
    
    cell_labels{n_dset} = max_idx(sel_idx);
    data_all{n_dset} = cat(2,resp_vals{:});
end

data_all2 = cat(1,data_all{:});
hasnan1 = logical(sum(isnan(data_all2),2));
data_all2 = data_all2(~hasnan1,:);

cell_labels2 = cat(1,cell_labels{:});
cell_labels2 = cell_labels2(~hasnan1);

[lr_data, residual_var, residual_var_pca] = f_dv_run_dred2(data_all2', method, dist_metric);

num_dim = numel(lr_data);

title_tag2 = sprintf('%s; resp %s', title_tag1, resp_cell_sel);

figure; hold on
plot([0, 1:numel(residual_var_pca)], [1; residual_var_pca], 'o-', 'Linewidth', 2);
if ~strcmpi(method, 'pca')
    plot([0, 1:numel(residual_var)], [1; residual_var], 'o-', 'Linewidth', 2);
    legend('PCA', method)
end
title(sprintf('Residual variance proj cells; %s', title_tag2), 'interpreter', 'none');
xlabel('number components used');
ylabel('Residual variance');

if num_dim >= 2
    plot_data = lr_data{2};
    figure; hold on
    plot(plot_data(:,1), plot_data(:,2), 'ok', 'LineWidth', 1)
    for n_cell = 1:numel(cell_labels2)
        if cell_labels2(n_cell)
            color1 = app.ops.context_types_all_colors2{tn_all(cell_labels2(n_cell))};
        else
            color1 = [.6 .6 .6];
        end
        plot(plot_data(n_cell, 1), plot_data(n_cell, 2), '.', 'color', color1, 'LineWidth', 4, 'MarkerSize', 15)
    end
    title(sprintf('low rank proj cells 2d; %s',title_tag2), 'interpreter', 'none');
end

if num_dim >= 3
    plot_data = lr_data{3};
    figure; hold on
    plot3(plot_data(:,1), plot_data(:,2), plot_data(:,3), 'ok', 'LineWidth', 1)
    for n_cell = 1:numel(cell_labels2)
        if cell_labels2(n_cell)
            color1 = app.ops.context_types_all_colors2{tn_all(cell_labels2(n_cell))};
        else
            color1 = [.6 .6 .6];
        end
        %plot3([0 lr_data3d(n_tn, 1)], [0 lr_data3d(n_tn, 2)], [0 lr_data3d(n_tn, 3)], '-', 'color', app.ops.context_types_all_colors2{tn_all(n_tn)}, 'LineWidth', 2)
        plot3(plot_data(n_cell, 1), plot_data(n_cell, 2), plot_data(n_cell, 3), '.', 'color', color1, 'LineWidth', 4, 'MarkerSize', 15)
    end
    title(sprintf('low rank proj cells 3d; %s', title_tag2), 'interpreter', 'none');
    grid on;

    title_tag3 = sprintf('%s; resp %s; cell ave', title_tag2, resp_cell_sel);
    data_mean = zeros(numel(tn_all), 3);
    for n_nt = 1:numel(tn_all)
        idx1 = cell_labels2 == n_nt;
        data_mean(n_nt,:) = mean(lr_data{3}(idx1,:),1);
    end
    f_dv_plot3_pc2(data_mean, tn_all, [], title_tag3, app.ops.context_types_all_colors2)


    if strcmpi(method ,'PCA')
        if num_dim >= 6
            plot_data = lr_data{6};
            figure; hold on;
            plot3(plot_data(:,4), plot_data(:,5), plot_data(:,6), 'ok', 'LineWidth', 1)
            for n_cell = 1:numel(cell_labels2)
                if cell_labels2(n_cell)
                    color1 = app.ops.context_types_all_colors2{tn_all(cell_labels2(n_cell))};
                else
                    color1 = [.6 .6 .6];
                end
                %plot3([0 lr_data3d(n_tn, 1)], [0 lr_data3d(n_tn, 2)], [0 lr_data3d(n_tn, 3)], '-', 'color', app.ops.context_types_all_colors2{tn_all(n_tn)}, 'LineWidth', 2)
                plot3(plot_data(n_cell, 4), plot_data(n_cell, 5), plot_data(n_cell, 6), '.', 'color', color1, 'LineWidth', 4, 'MarkerSize', 15)
            end
            xlabel('PC4'); ylabel('PC5'); zlabel('PC6');
            title(sprintf('low rank proj cells 3d; %s', title_tag2), 'interpreter', 'none');
            grid on;

            data_mean = zeros(numel(tn_all), 3);
            for n_nt = 1:numel(tn_all)
                idx1 = cell_labels2 == tn_all(n_nt);
                data_mean(n_nt,:) = mean(lr_data{6}(idx1,4:6),1);
            end
            f_dv_plot3_pc2(data_mean, tn_all, [], title_tag3, app.ops.context_types_all_colors2)
            xlabel('PC4'); ylabel('PC5'); zlabel('PC6');
        end
    end
end

end