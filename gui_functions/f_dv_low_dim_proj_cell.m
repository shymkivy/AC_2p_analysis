function f_dv_low_dim_proj_cell(app)

method = app.DimredmethodDropDown.Value;
dist_metric = app.DistmethodDropDown.Value; % pca isomap

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

shadow_axis_locs = [app.ReflectXCheckBox.Value, app.ReflectYCheckBox.Value, app.ReflectZCheckBox.Value] + 1;

num_dsets = size(data,1);

[tn0, con_idx] = f_dv_get_trial_number(app, data.MMN_freq{1});
tn00 = tn0(1,:);

[num_gr, num_tn] = size(tn0);

title_tag1 = sprintf('%s; %s', title_tag, method);
if strcmpi(method, 'isomap')
    title_tag1 = sprintf('%s; dist %s', title_tag1, dist_metric);
end

if strcmpi(app.ResponsivecellsselectDropDown.Value, 'All')
    resp_cell_sel = 'All';
else
    resp_cell_sel = 'Resp marg';
end

cell_labels = cell(num_dsets,num_gr);
data_all = cell(num_dsets,num_gr);
for n_dset = 1:num_dsets
    ddata = data(n_dset,:);
    stats1 = cat(1,ddata.stats{n_pl});
    tn1 = f_dv_get_trial_number(app, ddata.MMN_freq{1});

    for n_gr = 1:num_gr
        [sel_cells, resp_vals, resp_vals2, ~, resp_cells] = f_dv_get_resp_vals_cells(app, stats1, tn1(n_gr,:), [], resp_cell_sel);
        
        sel_idx = logical(sum(sel_cells,2));
        resp_idx = logical(sum(resp_cells,2));
        [~, max_idx] = max(resp_vals2, [], 2);
        
        max_idx(~resp_idx) = 0;

        cell_labels{n_dset, n_gr} = max_idx(sel_idx);
        data_all{n_dset, n_gr} = cat(2,resp_vals{:});
    end
end

data_all2 = cat(1,data_all{:});
[data_all2, rem_idx] = f_dv_fix_nan_trials(data_all2, app.nanhandlemetDropDown.Value);

cell_labels2 = cat(1,cell_labels{:});
cell_labels2 = cell_labels2(~rem_idx);


params.method = app.DimredmethodDropDown.Value;
params.dist_metric = app.DistmethodDropDown.Value;
params.subtract_mean = app.subtractmeanCheckBox.Value;
params.scale_by_var = app.scalebyvarCheckBox.Value;
params.plot_subtrat_mean = app.plotsubmeanCheckBox.Value;
[lr_data, residual_var, residual_var_pca, subtr_mean] = f_dv_run_dred(data_all2', params);

num_dim = size(lr_data,2);

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
    lr_data2d = lr_data(:,1:2);
    figure; hold on
    plot(lr_data2d(:,1), lr_data2d(:,2), 'ok', 'LineWidth', 1)
    for n_cell = 1:numel(cell_labels2)
        if cell_labels2(n_cell)
            color1 = app.ops.context_types_all_colors2{tn00(cell_labels2(n_cell))};
        else
            color1 = [.6 .6 .6];
        end
        plot(lr_data2d(n_cell, 1), lr_data2d(n_cell, 2), '.', 'color', color1, 'LineWidth', 4, 'MarkerSize', 15)
    end
    title(sprintf('low rank proj cells 2d; %s',title_tag2), 'interpreter', 'none');
end

if num_dim >= 3
    lr_data3d = lr_data(:,1:3);
    figure; hold on
    plot3(lr_data3d(:,1), lr_data3d(:,2), lr_data3d(:,3), 'ok', 'LineWidth', 1)
    for n_cell = 1:numel(cell_labels2)
        if cell_labels2(n_cell)
            color1 = app.ops.context_types_all_colors2{tn00(cell_labels2(n_cell))};
        else
            color1 = [.6 .6 .6];
        end
        %plot3([0 lr_data3d(n_tn, 1)], [0 lr_data3d(n_tn, 2)], [0 lr_data3d(n_tn, 3)], '-', 'color', app.ops.context_types_all_colors2{tn_all(n_tn)}, 'LineWidth', 2)
        plot3(lr_data3d(n_cell, 1), lr_data3d(n_cell, 2), lr_data3d(n_cell, 3), '.', 'color', color1, 'LineWidth', 4, 'MarkerSize', 15)
    end
    title(sprintf('low rank proj cells 3d; %s', title_tag2), 'interpreter', 'none');
    grid on;

    title_tag3 = sprintf('%s; resp %s; cell ave', title_tag2, resp_cell_sel);
    data_mean = zeros(num_tn, 3);
    for n_nt = 1:num_tn
        idx1 = cell_labels2 == n_nt;
        data_mean(n_nt,:) = mean(lr_data3d(idx1,:),1);
    end
    f_dv_plot3_pc(data_mean, tn00, [], title_tag3, app.ops.context_types_all_colors2, [], con_idx, shadow_axis_locs)


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
            f_dv_plot3_pc(data_mean, tn_all, [], title_tag3, app.ops.context_types_all_colors2)
            xlabel('PC4'); ylabel('PC5'); zlabel('PC6');
        end
    end
end

end