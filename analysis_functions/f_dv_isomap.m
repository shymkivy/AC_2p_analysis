function f_dv_isomap(app)

params = f_dv_gather_params(app);
ops = app.ops;
ddata = app.ddata;
cdata = f_dv_get_cdata(app);

volume_period = ddata.proc_data{1}.frame_data.volume_period;
more_smooth_SD = 0;

n_pl = params.n_pl;

tn_all = f_dv_get_trial_number(params);
tt_all = ops.context_types_all(tn_all)';

stim_times = ddata.stim_frame_index{n_pl};
mmn_freq = ddata.MMN_freq{1};

trial_types = ddata.trial_types{1};

[~, trial_frames] = f_dv_compute_window_t(params.trial_window, cdata.volume_period);

firing_rate = cat(1,cdata.S_sm);
num_cells = size(firing_rate,1);

if more_smooth_SD
    firing_rate = f_smooth_gauss(firing_rate, more_smooth_SD/volume_period);
end
if params.shuffle_cells
    firing_rate = firing_rate(randperm(num_cells),:);
end

trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trial_frames);
[trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, ops);

num_t = size(trial_data_sort,2);

trial_idx = logical(sum(tt_all == trial_types_wctx,2));
trial_seq = trial_types_wctx(trial_idx);
trial_data_sort2 = trial_data_sort_wctx(:,:,trial_idx);

if params.shuffle_trials
    num_trials = sum(trial_idx);
    shuff_seq = randperm(num_trials);
    trial_data_sort2 = trial_data_sort2(:,:,shuff_seq);
    trial_seq = trial_seq(shuff_seq);
end

tn_seq = f_tt_to_tn(trial_seq, ops,0);


if params.sort_by_trial_type
    [tn_seq, sort_idx] = sort(tn_seq);
    trial_data_sort2 = trial_data_sort2(:,:,sort_idx);
    trial_seq = trial_seq(sort_idx);
end

firing_rate2 = reshape(trial_data_sort2, num_cells, []);

% remove inactive cells
active_cells = sum(firing_rate2,2) ~= 0;
firing_rate2(~active_cells,:) = [];

if params.sort_by_similarity
    hc_params.method = 'ward'; % ward(inner square), average, single(shortest)
    hc_params.distance_metric = 'cosine'; % none, euclidean, squaredeuclidean, cosine, hammilarity
    hc_params.plot_dist_mat = params.plot_sorting_stuff;
    hc_params.plot_clusters = params.plot_sorting_stuff;
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