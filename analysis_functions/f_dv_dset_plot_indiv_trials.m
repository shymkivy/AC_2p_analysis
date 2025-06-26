function f_dv_dset_plot_indiv_trials(app)
%% plot raster of firing rates without enseble analysis

% app stuff
ops = app.ops;
params = f_dv_gather_params(app);

tn_all = f_dv_get_trial_number(params);
ddata = app.ddata;
cdata = f_dv_get_cdata(app);

method = params.dim_red_method;
dist_metric = params.distance_method; % pca isomap

mean_win = [.1 .9];

tt_all = ops.context_types_all(tn_all)';

firing_rate = cat(1,cdata.S_sm);
stats1 = cat(1,ddata.stats{:});

title_tag1 = sprintf('%s; %s', ddata.dset_name_full{1}, method);
if strcmpi(method, 'isomap')
    title_tag1 = sprintf('%s; dist %s', title_tag1, dist_metric);
end

[sel_cells] = f_dv_get_resp_vals_cells(stats1, tn_all, params);
sel_cells2 = logical(sum(sel_cells,2));

firing_rate = firing_rate(sel_cells2,:);

num_cells = sum(sel_cells2);

if isempty(tn_all)
    disp('Analyzing full trace')
    firing_rate3 = firing_rate;
    tn_seq_plot = [];
else
    if isempty(tt_all)
        disp('Specified trial type does not exist');
    else
        stim_times = ddata.stim_frame_index{1};
        mmn_freq = ddata.MMN_freq{1};
        
        % params.trial_window
        [trial_t, trial_frames] = f_dv_compute_window_t(mean_win, cdata(1).volume_period);

        trial_types = ddata.trial_types{1};
        trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trial_frames);
         
        if ~isempty(mmn_freq)
            [trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, ops);
        else
            trial_data_sort_wctx = trial_data_sort;
            trial_types_wctx = trial_types;
        end
        
        trial_data_sort_wctx_2d = reshape(mean(trial_data_sort_wctx,2), num_cells, []);
        
        trial_idx = logical(sum(tt_all == trial_types_wctx,2));
        trial_seq = trial_types_wctx(trial_idx);
        trial_data_sort2 = trial_data_sort_wctx_2d(:,trial_idx);

        tn_seq = f_tt_to_tn(trial_seq, ops,0);

%         if app.sortbytrialtypeCheckBox.Value
%             [tn_seq, sort_idx] = sort(tn_seq);
%             trial_data_sort2 = trial_data_sort2(:,:,sort_idx);
%             trial_seq = trial_seq(sort_idx);
%         end
    end
end

%% clustering and plotting

[lr_data, residual_var, residual_var_pca] = f_dv_run_dred2(trial_data_sort2, method, dist_metric);

title_tag2 = sprintf('%s; resp %s', title_tag1, params.responsive_cells_select);

figure; hold on
plot([0, 1:numel(residual_var_pca)], [1; residual_var_pca], 'o-', 'Linewidth', 2);
if ~strcmpi(method, 'pca')
    plot([0, 1:numel(residual_var)], [1; residual_var], 'o-', 'Linewidth', 2);
    legend('PCA', method)
end
title(sprintf('Residual variance proj cells; %s', title_tag2), 'interpreter', 'none');
xlabel('number components used');
ylabel('Residual variance');

data_plot = lr_data{2};
figure; hold on
plot(data_plot(:,1), data_plot(:,2), 'ok', 'Linewidth', 1)
for n_tn = 1:numel(tn_seq)
    color1 = ops.context_types_all_colors2{tn_seq(n_tn)};
    plot(data_plot(n_tn, 1), data_plot(n_tn, 2), '.', 'color', color1, 'LineWidth', 2, 'MarkerSize', 20)
end
title(sprintf('low rank proj trials 2d; %s', title_tag2), 'interpreter', 'none');

data_plot = lr_data{3};
figure; hold on
plot3(data_plot(:,1), data_plot(:,2), data_plot(:,3), 'ok', 'LineWidth', 1)
for n_tn = 1:numel(tn_seq)
    color1 = ops.context_types_all_colors2{tn_seq(n_tn)};
    %plot3([0 lr_data3d(n_tn, 1)], [0 lr_data3d(n_tn, 2)], [0 lr_data3d(n_tn, 3)], '-', 'color', app.ops.context_types_all_colors2{tn_all(n_tn)}, 'LineWidth', 2)
    plot3(data_plot(n_tn, 1), data_plot(n_tn, 2), data_plot(n_tn, 3), '.', 'color', color1, 'LineWidth', 4, 'MarkerSize', 15)
end
title(sprintf('low rank proj cells 3d; %s', title_tag2), 'interpreter', 'none');
grid on;

end