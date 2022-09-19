function f_dv_plot_ens_raster(app)
%% needs ens data

n_tr = 6;

tn_all = f_dv_get_trial_number(app);
tt_all = app.ops.context_types_all(tn_all)';

ddata = app.ddata;
cdata = f_dv_get_cdata(app);
firing_rate = cat(1,cdata.S_sm);

ens_stats = ddata.ensemble_stats{1};
ens_tuning = ddata.ensemble_tuning_stats{1};
ens_data = ddata.ensembles{1};

raster_lr = ens_data.coeffs*ens_data.scores;
figure; imagesc(raster_lr)

figure; imagesc(firing_rate)

accepted_ens = ens_stats.accepted_ensembles;

ens_cells_list = ens_data.cells.ens_list(accepted_ens);
ens_cell_coeffs = ens_data.coeffs(:,accepted_ens);
scores = ens_data.scores(accepted_ens,:);

resp_ens = find(ens_tuning.resp_cells_peak(:,n_tr));

figure; hold on;
for n_ens = 1:numel(resp_ens)
   n_ens2 =  resp_ens(n_ens);
   plot(scores(n_ens2,:))
end

cell_ens_idx = false(num_cells, numel(resp_ens));
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


stim_frame_index = ddata.stim_frame_index{1}(ddata.trial_types{1} == n_tr);

trial_num_baseline_resp_frames = ddata.trial_window{1}.trial_num_baseline_resp_frames;
trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_frame_index, trial_num_baseline_resp_frames);

trial_data_sort_2d = reshape(trial_data_sort, num_cells, []);
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