function f_dv_plot_raster(app)

n_pl = app.mplSpinner.Value;
tn_all = f_dv_get_trial_number(app);
tt_all = app.ops.context_types_all(tn_all)';

stim_times = app.ddata.stim_frame_index{n_pl};
mmn_freq = app.ddata.MMN_freq{1};
trig_window = app.working_ops.trial_num_baseline_resp_frames;
trial_types = app.ddata.trial_types{1};

cdata = f_dv_get_cdata(app);

num_cells = sum([cdata.num_cells]);
firing_rate = cat(1,cdata.S);

if app.shufflecellsCheckBox.Value
    firing_rate = firing_rate(randperm(num_cells),:);
end

trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trig_window);

if ~isempty(mmn_freq)
    [trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, app.ops);
else
    trial_data_sort_wctx = trial_data_sort;
    trial_types_wctx = trial_types;
end

num_t = size(trial_data_sort,2);

trial_idx = logical(sum(tt_all == trial_types_wctx,2));
trial_seq = trial_types_wctx(trial_idx);
trial_data_sort2 = trial_data_sort_wctx(:,:,trial_idx);

if app.shuffletrialsCheckBox.Value
    num_trials = sum(trial_idx);
    shuff_seq = randperm(num_trials);
    trial_data_sort2 = trial_data_sort2(:,:,shuff_seq);
    trial_seq = trial_seq(shuff_seq);
end

tn_seq = f_tt_to_tn(trial_seq, app.ops,0);


if app.sortbytrialtypeCheckBox.Value
    [tn_seq, sort_idx] = sort(tn_seq);
    trial_data_sort2 = trial_data_sort2(:,:,sort_idx);
    trial_seq = trial_seq(sort_idx);
end

firing_rate2 = reshape(trial_data_sort2, num_cells, []);

% remove inactive cells
active_cells = sum(firing_rate2,2) ~= 0;
firing_rate2(~active_cells,:) = [];

if app.sortbycellsimilarityCheckBox.Value
    hc_params.method = 'ward'; % ward(inner square), average, single(shortest)
    hc_params.distance_metric = 'cosine'; % none, euclidean, squaredeuclidean, cosine, hammilarity
    hc_params.plot_dist_mat = app.plotsortingstuffCheckBox.Value;
    hc_params.plot_clusters = app.plotsortingstuffCheckBox.Value;
    hc_params.num_clust = [];
    hc_params.title_tag = 'Coeffs (cells)';
    hclust_out_cell = f_hcluster_wrap(firing_rate2, hc_params);
    ord_cell = hclust_out_cell.dend_order;
else
    ord_cell = 1:num_cells;
end

tn_seq_plot = reshape(repmat(tn_seq, [1 num_t])',[],1);
f_plot_raster_mean(firing_rate2(ord_cell,:), 1, tn_seq_plot, app.ops.context_types_all_colors2)

end