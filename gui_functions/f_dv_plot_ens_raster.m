function f_dv_plot_ens_raster(app)

plot_selected_tr = 1;

n_tr = f_dv_get_trial_number(app);
firing_rate = app.cdata.S;

ens_stats = app.ddata.ensemble_stats{1};
ens_tuning = app.ddata.ensemble_tuning{1};
ens_data = app.ddata.ensembles{1}.ens_out;

%raster_lr = ens_data.coeffs*ens_data.scores;
%figure; imagesc(raster_lr)
%figure; imagesc(firing_rate)

if isfield(ens_stats, 'accepted_ensembles')
    accepted_ens = ens_stats.accepted_ensembles;
    ens_cells_list = ens_data.cells.ens_list(accepted_ens);
    ens_cell_coeffs = ens_data.coeffs(:,accepted_ens);
    scores = ens_data.scores(accepted_ens,:);
else
    ens_cells_list = ens_data.cells.ens_list;
    ens_cell_coeffs = ens_data.coeffs;
    scores = ens_data.scores;
end



resp_ens = find(logical(sum(ens_tuning.cell_is_resp(:,n_tr),2)));
% 
% figure; hold on;
% for n_ens = 1:numel(resp_ens)
%    n_ens2 =  resp_ens(n_ens);
%    plot(scores(n_ens2,:))
% end

cell_ens_idx = false(app.cdata.num_cells, numel(resp_ens));
for n_ens = 1:numel(resp_ens)
    n_ens2 =  resp_ens(n_ens);
    cell_list = ens_cells_list{n_ens2};
    cell_ens_idx(cell_list,n_ens) = 1;
end

% for n_ens = 1:numel(resp_ens)
%     n_ens2 =  resp_ens(n_ens);
%     cell_list = ens_cells_list{n_ens2};
%     [~, sort_idx] = sort(ens_cell_coeffs(cell_list,n_ens2), 'descend');
%     figure;
%     ax1 = subplot(2,1,1);
%     imagesc(firing_rate(cell_list(sort_idx),:));
%     title(['ens ' num2str(n_ens2)]);
%     ax2 = subplot(2,1,2); hold on;
%     plot(sum(firing_rate(cell_list(sort_idx),:))/ max(sum(firing_rate(cell_list(sort_idx),:))));
%     plot(scores(n_ens2,:)/max(scores(n_ens2,:)))
%     linkaxes([ax1,ax2],'x');
% end

all_cell_list = [];
for n_ens = 1:numel(resp_ens)
    n_ens2 =  resp_ens(n_ens);
    cell_list = ens_cells_list{n_ens2};
    [~, sort_idx] = sort(ens_cell_coeffs(cell_list,n_ens2), 'descend');
    all_cell_list = [all_cell_list; cell_list(sort_idx)];
end

all_cell_list_uniq = unique(all_cell_list, 'stable');
num_gr_cells = numel(all_cell_list_uniq);
firing_rate_ens = firing_rate(all_cell_list_uniq,:);

if app.SortenstrialsCheckBox.Value
    trial_types = app.ddata.trial_types{1};
    stim_frame_index = app.ddata.stim_frame_index{1};
    trial_num_baseline_resp_frames = app.ddata.trial_window{1}.trial_num_baseline_resp_frames;
    trial_data_sort_ens = f_get_stim_trig_resp(firing_rate_ens, stim_frame_index, trial_num_baseline_resp_frames);
    trial_data_sort_scores = f_get_stim_trig_resp(scores(resp_ens,:), stim_frame_index, trial_num_baseline_resp_frames);
    
    
    if plot_selected_tr
        [trial_data_sort_ens_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort_ens, trial_types, app.ddata.MMN_freq{1}, app.ops);
        sel_tr = logical(sum(trial_types_wctx==app.ops.context_types_all(n_tr)',2));
        trial_data_sort_ens2 = trial_data_sort_ens_wctx(:,:,sel_tr);
        trial_types2 = trial_types_wctx(sel_tr);
        
        [trial_data_sort_scores_wctx, ~] =  f_s3_add_ctx_trials(trial_data_sort_scores, trial_types, app.ddata.MMN_freq{1}, app.ops);
        trial_data_sort_scores2 = trial_data_sort_scores_wctx(:,:,sel_tr);
    else
        trial_data_sort_ens2 = trial_data_sort_ens;
        trial_data_sort_scores2 = trial_data_sort_scores;
        trial_types2 = trial_types;
    end
    [~, idx1] = sort(trial_types2);
    
    trial_data_sort_ens4 = trial_data_sort_ens2(:,:,idx1);
    
    num_tr = size(trial_data_sort_scores2,3);
    trial_data_sort_ens5 = reshape(trial_data_sort_scores2(:,:,idx1), [], num_tr);
    params.plot_stuff = 0;
    hclust_out = f_hcluster_wrap(trial_data_sort_ens5', params);
    
    trial_data_sort_ens5 = trial_data_sort_ens4(:,:,hclust_out.dend_order);
    
    trial_data_sort_ens3 = reshape(trial_data_sort_ens5, num_gr_cells, []);
else
    trial_data_sort_ens3 = firing_rate_ens;
end


figure; 
subplot(1,10,1:9);
imagesc(trial_data_sort_ens3);
subplot(1,10,10);
imagesc(cell_ens_idx(all_cell_list_uniq,:));

end