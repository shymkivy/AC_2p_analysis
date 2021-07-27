function f_dv_plot_ens_deets(app)

%%
n_tr = 30;
firing_rate = app.cdata.S;

ens_stats = app.ddata.ensemble_stats{1};

ens_tuning = app.ddata.ensemble_tuning{1};

ens_data = app.ddata.ensembles{1}.ens_out;

%raster_lr = ens_data.coeffs*ens_data.scores;
%figure; imagesc(raster_lr)
%figure; imagesc(firing_rate)

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


trial_types = app.ddata.trial_types{1};
stim_frame_index = app.ddata.stim_frame_index{1};
trial_num_baseline_resp_frames = app.ddata.trial_window{1}.trial_num_baseline_resp_frames;

trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_frame_index, trial_num_baseline_resp_frames);
[trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, app.ddata.MMN_freq{1}, app.ops);


trial_data_sort2 = trial_data_sort_wctx(:,:,trial_types_wctx==app.ops.context_types_all(n_tr));

num_gr_cells = numel(all_cell_list_uniq)
peak_vals = squeeze(max(trial_data_sort2(all_cell_list_uniq,:,:), [], 2));

params.plot_stuff = 1;
hclust_out = f_hcluster_wrap(peak_vals', params);

trial_data_sort_2d = reshape(trial_data_sort2(:,:,hclust_out.dend_order), app.cdata.num_cells, []);
raster_ens = trial_data_sort_2d(all_cell_list_uniq,:);

figure; 
subplot(1,10,1:9);
imagesc(raster_ens);
subplot(1,10,10);
imagesc(cell_ens_idx(all_cell_list_uniq,:));

[~, idx1] = sort(trial_types);
trial_data_sort_sort = trial_data_sort(:,:,idx1);
trial_data_sort_sort_2d = reshape(trial_data_sort_sort, app.cdata.num_cells, []);
raster_ens_all_tr = trial_data_sort_sort_2d(all_cell_list_uniq,:);

figure; 
subplot(1,10,1:9);
imagesc(raster_ens_all_tr);
subplot(1,10,10);
imagesc(cell_ens_idx(all_cell_list_uniq,:));



D = f_pdist2_YS(raster_ens, raster_ens, 'cosine');

figure; imagesc(1-D)

figure; imagesc(scores(resp_ens,:))

figure; imagesc(scores(resp_ens,:))

figure; imagesc(scores)

trial_scores_sort = f_get_stim_trig_resp(scores, stim_frame_index, trial_num_baseline_resp_frames);

trial_scores_sort_2d = reshape(trial_scores_sort(:,:,idx1), size(scores,1), []);
figure; imagesc(trial_scores_sort_2d)


end