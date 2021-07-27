function f_dv_plot_raster(app)

n_pl = app.mplSpinner.Value;
tn_all = f_dv_get_trial_number(app);
tt_all = app.ops.context_types_all(tn_all)';

stim_times = app.ddata.stim_frame_index{n_pl};
mmn_freq = app.ddata.MMN_freq{1};
trig_window = app.working_ops.trial_num_baseline_resp_frames;
trial_types = app.ddata.trial_types{1};

firing_rate = app.cdata.S;
num_cells = size(firing_rate,1);

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

tn_seq_plot = reshape(repmat(tn_seq, [1 num_t])',[],1);
f_plot_raster_mean(firing_rate2(ord_cell,:), 1, tn_seq_plot, app.ops.context_types_all_colors2)

%%
n_tr = 6;
firing_rate = app.cdata.S;

ens_stats = app.ddata.ensemble_stats{1};

ens_tuning = app.ddata.ensemble_tuning{1};

ens_data = app.ddata.ensembles{1}.ens_out;

raster_lr = ens_data.coeffs*ens_data.scores;
figure; imagesc(raster_lr)

figure; imagesc(firing_rate)

accepted_ens = ens_stats.accepted_ensembles;

ens_cells_list = ens_data.cells.ens_list(accepted_ens);
ens_cell_coeffs = ens_data.coeffs(:,accepted_ens);
scores = ens_data.scores(accepted_ens,:);

resp_ens = find(ens_tuning.cell_is_resp(:,n_tr));

figure; hold on;
for n_ens = 1:numel(resp_ens)
   n_ens2 =  resp_ens(n_ens);
   plot(scores(n_ens2,:))
end

cell_ens_idx = false(app.cdata.num_cells, numel(resp_ens));
for n_ens = 1:numel(resp_ens)
    n_ens2 =  resp_ens(n_ens);
    cell_list = ens_cells_list{n_ens2};
    cell_ens_idx(cell_list,n_ens) = 1;
end

for n_ens = 1:numel(resp_ens)
    n_ens2 =  resp_ens(n_ens);
    cell_list = ens_cells_list{n_ens2};
    [~, sort_idx] = sort(ens_cell_coeffs(cell_list,n_ens2), 'descend');
    figure;
    ax1 = subplot(2,1,1);
    imagesc(firing_rate(cell_list(sort_idx),:));
    title(['ens ' num2str(n_ens2)]);
    ax2 = subplot(2,1,2); hold on;
    plot(sum(firing_rate(cell_list(sort_idx),:))/ max(sum(firing_rate(cell_list(sort_idx),:))));
    plot(scores(n_ens2,:)/max(scores(n_ens2,:)))
    linkaxes([ax1,ax2],'x');
end

all_cell_list = [];
for n_ens = 1:numel(resp_ens)
    n_ens2 =  resp_ens(n_ens);
    cell_list = ens_cells_list{n_ens2};
    [~, sort_idx] = sort(ens_cell_coeffs(cell_list,n_ens2), 'descend');
    all_cell_list = [all_cell_list; cell_list(sort_idx)];
end

all_cell_list_uniq = unique(all_cell_list, 'stable');

figure; 
subplot(1,10,1:9);
imagesc(firing_rate(all_cell_list_uniq,:));
subplot(1,10,10);
imagesc(cell_ens_idx(all_cell_list_uniq,:));


stim_frame_index = app.ddata.stim_frame_index{1}(app.ddata.trial_types{1} == n_tr);

trial_num_baseline_resp_frames = app.ddata.trial_window{1}.trial_num_baseline_resp_frames;
trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_frame_index, trial_num_baseline_resp_frames);

trial_data_sort_2d = reshape(trial_data_sort, app.cdata.num_cells, []);
raster_ens = trial_data_sort_2d(all_cell_list_uniq,:);

figure; 
subplot(1,10,1:9);
imagesc(raster_ens);
subplot(1,10,10);
imagesc(cell_ens_idx(all_cell_list_uniq,:));

figure; imagesc(raster_ens)

D = f_pdist2_YS(raster_ens, raster_ens, 'cosine');

figure; imagesc(1-D)




end