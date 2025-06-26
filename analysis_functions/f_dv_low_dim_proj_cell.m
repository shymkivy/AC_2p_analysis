function f_dv_low_dim_proj_cell(app)

% app stuff
ops = app.ops;
params = f_dv_gather_params(app);

[data, title_tag] = f_dv_get_data_by_mouse_selection(app.data, params);
%[region_num, reg_tag, leg_list] = f_dv_get_region_sel_val(params, ops);
[tn0, con_idx] =  f_dv_get_trial_number(params, data.MMN_freq{1});

num_dsets = size(data,1);

tn00 = tn0(1,:);

[num_gr, num_tn] = size(tn0);

title_tag1 = sprintf('%s; %s', title_tag, params.dim_red_method);
if strcmpi(params.dim_red_method, 'isomap')
    title_tag1 = sprintf('%s; dist %s', title_tag1, params.distance_method);
end

if ~strcmpi(params.responsive_cells_select, 'All')
    params.responsive_cells_select = 'Resp marg';
end

cell_labels = cell(num_dsets,num_gr);
data_all = cell(num_dsets,num_gr);
for n_dset = 1:num_dsets
    ddata = data(n_dset,:);
    stats1 = cat(1,ddata.stats{params.planes});

    tn1 = f_dv_get_trial_number(params, ddata.MMN_freq{1});

    for n_gr = 1:num_gr
        [sel_cells, resp_vals, resp_vals2, ~, resp_cells] = f_dv_get_resp_vals_cells(stats1, tn1(n_gr,:), params);
        
        sel_idx = logical(sum(sel_cells,2));
        resp_idx = logical(sum(resp_cells,2));
        [~, max_idx] = max(resp_vals2, [], 2);
        
        max_idx(~resp_idx) = 0;

        cell_labels{n_dset, n_gr} = max_idx(sel_idx);
        data_all{n_dset, n_gr} = cat(2,resp_vals{:});
    end
end

data_all2 = cat(1,data_all{:});
[data_all2, rem_idx] = f_dv_fix_nan_trials(data_all2, params.non_handle_method);

cell_labels2 = cat(1,cell_labels{:});
cell_labels2 = cell_labels2(~rem_idx);

dred_params.method = params.dim_red_method;
dred_params.dist_metric = params.distance_method;
dred_params.subtract_mean = params.subtract_mean;
dred_params.scale_by_var = params.scale_by_var;
dred_params.plot_subtrat_mean = params.plot_subtracted_mean;
[lr_data, residual_var, residual_var_pca, ~] = f_dv_run_dred(data_all2', dred_params);

num_dim = size(lr_data,2);

title_tag2 = sprintf('%s; resp %s', title_tag1, params.responsive_cells_select);

figure; hold on
plot([0, 1:numel(residual_var_pca)], [1; residual_var_pca], 'o-', 'Linewidth', 2);
if ~strcmpi(params.dim_red_method, 'pca')
    plot([0, 1:numel(residual_var)], [1; residual_var], 'o-', 'Linewidth', 2);
    legend('PCA', params.dim_red_method)
end
title(sprintf('Residual variance proj cells; %s', title_tag2), 'interpreter', 'none');
xlabel('number components used');
ylabel('Residual variance');

if strcmpi(params.num_axes_plot, '2')
    num_plots = ceil(params.num_comp_plot/2);
    % 2 comp  bs
    for n_pl = 1:num_plots
        pcs = [(n_pl-1)*2+1, (n_pl-1)*2+2];
        title_tag3 = sprintf('low rank proj cells; pl%d; %s', n_pl, title_tag2);

        f_dv_plot2_pc_bytn(lr_data, tn00, cell_labels2, pcs, title_tag3, ops.context_types_all_colors2, params.render_painters)
    end
elseif strcmpi(params.num_axes_plot, '3')
    num_plots = ceil(params.num_comp_plot/3);
    % 3 comp  bs
    for n_pl = 1:num_plots
        pcs = [(n_pl-1)*3+1, (n_pl-1)*3+2, (n_pl-1)*3+3];
        title_tag3 = sprintf('low rank proj cells; pl%d; %s', n_pl, title_tag2);
        
        f_dv_plot3_pc_bytn(lr_data, tn00, cell_labels2, pcs, title_tag3, ops.context_types_all_colors2, params.shadow_on3d, params.shadow_axis_locs, params.render_painters, params.grid_on, params.reverse_xyz)
    end
end

% 
% if num_dim >= 2
%     lr_data2d = lr_data(:,1:2);
%     figure; hold on
%     plot(lr_data2d(:,1), lr_data2d(:,2), 'ok', 'LineWidth', 1)
%     for n_cell = 1:numel(cell_labels2)
%         if cell_labels2(n_cell)
%             color1 = app.ops.context_types_all_colors2{tn00(cell_labels2(n_cell))};
%         else
%             color1 = [.6 .6 .6];
%         end
%         plot(lr_data2d(n_cell, 1), lr_data2d(n_cell, 2), '.', 'color', color1, 'LineWidth', 4, 'MarkerSize', 15)
%     end
%     title(sprintf('low rank proj cells 2d; %s',title_tag2), 'interpreter', 'none');
% end
% 
% if num_dim >= 3
%     lr_data3d = lr_data(:,1:3);
%     figure; hold on
%     plot3(lr_data3d(:,1), lr_data3d(:,2), lr_data3d(:,3), 'ok', 'LineWidth', 1)
%     idx1 = cell_labels2 == 0;
%     plot3(lr_data3d(idx1, 1), lr_data3d(idx1, 2), lr_data3d(idx1, 3), '.', 'color', app.ops.context_types_all_colors2{tn00(n_tn)}, 'LineWidth', 4, 'MarkerSize', 15);
%     for n_tn = 1:num_tn
%         idx1 = cell_labels2 == n_tn;
%         plot3(lr_data3d(idx1, 1), lr_data3d(idx1, 2), lr_data3d(idx1, 3), '.', 'color', app.ops.context_types_all_colors2{tn00(n_tn)}, 'LineWidth', 4, 'MarkerSize', 15);
%     end
% 
%     % for n_cell = 1:numel(cell_labels2)
%     %     if cell_labels2(n_cell)
%     %         color1 = app.ops.context_types_all_colors2{tn00(cell_labels2(n_cell))};
%     %     else
%     %         color1 = [.6 .6 .6];
%     %     end
%     %     %plot3([0 lr_data3d(n_tn, 1)], [0 lr_data3d(n_tn, 2)], [0 lr_data3d(n_tn, 3)], '-', 'color', app.ops.context_types_all_colors2{tn_all(n_tn)}, 'LineWidth', 2)
%     %     plot3(lr_data3d(n_cell, 1), lr_data3d(n_cell, 2), lr_data3d(n_cell, 3), '.', 'color', color1, 'LineWidth', 4, 'MarkerSize', 15)
%     % end
%     title(sprintf('low rank proj cells 3d; %s', title_tag2), 'interpreter', 'none');
%     grid on;
% 
%     title_tag3 = sprintf('%s; resp %s; cell ave', title_tag2, resp_cell_sel);
%     data_mean = zeros(num_tn, 3);
%     for n_nt = 1:num_tn
%         idx1 = cell_labels2 == n_nt;
%         data_mean(n_nt,:) = mean(lr_data3d(idx1,:),1);
%     end
%     f_dv_plot3_pc(data_mean, tn00, [], title_tag3, app.ops.context_types_all_colors2, con_idx, 1, shadow_axis_locs)
% 
% 
%     if strcmpi(method ,'PCA')
%         if num_dim >= 6
%             plot_data = lr_data{6};
%             figure; hold on;
%             plot3(plot_data(:,4), plot_data(:,5), plot_data(:,6), 'ok', 'LineWidth', 1)
%             for n_cell = 1:numel(cell_labels2)
%                 if cell_labels2(n_cell)
%                     color1 = app.ops.context_types_all_colors2{tn_all(cell_labels2(n_cell))};
%                 else
%                     color1 = [.6 .6 .6];
%                 end
%                 %plot3([0 lr_data3d(n_tn, 1)], [0 lr_data3d(n_tn, 2)], [0 lr_data3d(n_tn, 3)], '-', 'color', app.ops.context_types_all_colors2{tn_all(n_tn)}, 'LineWidth', 2)
%                 plot3(plot_data(n_cell, 4), plot_data(n_cell, 5), plot_data(n_cell, 6), '.', 'color', color1, 'LineWidth', 4, 'MarkerSize', 15)
%             end
%             xlabel('PC4'); ylabel('PC5'); zlabel('PC6');
%             title(sprintf('low rank proj cells 3d; %s', title_tag2), 'interpreter', 'none');
%             grid on;
% 
%             data_mean = zeros(numel(tn_all), 3);
%             for n_nt = 1:numel(tn_all)
%                 idx1 = cell_labels2 == tn_all(n_nt);
%                 data_mean(n_nt,:) = mean(lr_data{6}(idx1,4:6),1);
%             end
%             f_dv_plot3_pc(data_mean, tn_all, [], title_tag3, app.ops.context_types_all_colors2, 1, shadow_axis_locs)
%             xlabel('PC4'); ylabel('PC5'); zlabel('PC6');
%         end
%     end
% end

end