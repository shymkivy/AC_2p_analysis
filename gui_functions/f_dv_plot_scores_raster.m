function f_dv_plot_scores_raster(app)

plot_selected_tr = 1;

n_tr = f_dv_get_trial_number(app);
ens_stats = app.ddata.ensemble_stats{1};
ens_tuning = app.ddata.ensemble_tuning{1};
ens_data = app.ddata.ensembles{1}.ens_out;

%raster_lr = ens_data.coeffs*ens_data.scores;
%figure; imagesc(raster_lr)
%figure; imagesc(firing_rate)

if isfield(ens_stats, 'accepted_ensembles')
    accepted_ens = ens_stats.accepted_ensembles;
    scores = ens_data.scores(accepted_ens,:);
else
    scores = ens_data.scores;
end

resp_ens = find(logical(sum(ens_tuning.cell_is_resp(:,n_tr),2)));


num_gr_scores = numel(resp_ens);
firing_rate_ens = scores(resp_ens,:);

if app.SortenstrialsCheckBox.Value
    trial_types = app.ddata.trial_types{1};
    stim_frame_index = app.ddata.stim_frame_index{1};
    trial_num_baseline_resp_frames = app.ddata.trial_window{1}.trial_num_baseline_resp_frames;
    trial_data_sort_ens = f_get_stim_trig_resp(firing_rate_ens, stim_frame_index, trial_num_baseline_resp_frames);
    
    if plot_selected_tr
        [trial_data_sort_ens_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort_ens, trial_types, app.ddata.MMN_freq{1}, app.ops);
        sel_tr = logical(sum(trial_types_wctx==app.ops.context_types_all(n_tr)',2));
        trial_data_sort_ens2 = trial_data_sort_ens_wctx(:,:,sel_tr);
        trial_types2 = trial_types_wctx(sel_tr);
    else
        trial_data_sort_ens2 = trial_data_sort_ens;
        trial_types2 = trial_types;
    end
    [~, idx1] = sort(trial_types2);
    trial_data_sort_ens4 = trial_data_sort_ens2(:,:,idx1);
    
    num_tr = size(trial_data_sort_ens4,3);
    trial_data_sort_ens5 = reshape(trial_data_sort_ens4, [], num_tr);
    params.plot_stuff = 0;
    hclust_out = f_hcluster_wrap(trial_data_sort_ens5', params);
    trial_data_sort_ens5 = trial_data_sort_ens4(:,:,hclust_out.dend_order);
    
    trial_data_sort_ens3 = reshape(trial_data_sort_ens5, num_gr_scores, []);
else
    trial_data_sort_ens3 = firing_rate_ens;
end

figure; 
imagesc(trial_data_sort_ens3);




end