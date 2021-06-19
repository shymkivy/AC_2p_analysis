function f_dv_isomap(app)

volume_period = app.ddata.proc_data{1}.frame_data.volume_period;
smooth_SD = 100;

n_pl = app.mplSpinner.Value;
tn_all = f_dv_get_trial_number(app);
tt_all = app.ops.context_types_all(tn_all)';

stim_times = app.ddata.stim_frame_index{n_pl};
mmn_freq = app.ddata.MMN_freq{1};
trig_window = app.working_ops.trial_num_baseline_resp_frames;
trial_types = app.ddata.trial_types{1};

firing_rate = app.cdata.S;
num_cells = size(firing_rate,1);

if smooth_SD
    firing_rate = f_smooth_gauss(firing_rate, smooth_SD/volume_period);
end
if app.shufflecellsCheckBox.Value
    firing_rate = firing_rate(randperm(num_cells),:);
end

trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trig_window);
[trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, app.ops);

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

no_dims = round(intrinsic_dim(firing_rate2, 'GMST'));


[X, labels] = generate_data('helix', 2000);
figure, scatter3(X(:,1), X(:,2), X(:,3), 5, labels); title('Original dataset'), drawnow
no_dims = round(intrinsic_dim(X, 'MLE'));
disp(['MLE estimate of intrinsic dimensionality: ' num2str(no_dims)]);

no_dims = 2

[mappedX, mapping] = compute_mapping(firing_rate2, 'PCA', 3);	
figure, scatter3(mappedX(:,1), mappedX(:,2), mappedX(:,3), 5); title('Result of PCA');

[mappedX, mapping] = compute_mapping(firing_rate2', 'Isomap', 2);	
figure, scatter(mappedX(:,1), mappedX(:,2), 5); title('Result of Isomap');

[mappedX, mapping] = compute_mapping(firing_rate2, 'GPLVM', no_dims, 7);	
figure, scatter(mappedX(:,1), mappedX(:,2), 5, labels); title('Result of GPLVM'); drawnow % (mapping.conn_comp)


end