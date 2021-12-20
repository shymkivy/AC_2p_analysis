function data_out = f_dv_ensless_dim_est(app, params)
% this will compute dimensionality of trials
params2 = f_dv_ensemble_params([]);
est_params_pca = params2.est_params_pca;
%%

firing_rate = cat(1,params.cdata.S_sm);

trial_types = params.ddata.trial_types{1};
stim_times = params.ddata.stim_frame_index{1};
trig_window = app.working_ops.trial_num_baseline_resp_frames;
mmn_freq = params.ddata.MMN_freq{1};

trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trig_window);
[trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, app.ops);

%%

stats1 = params.ddata.stats{1};

select_resp_cells = 1;

dim_est_tr = zeros(numel(app.ops.context_types_all),1);
num_cells_all = zeros(numel(app.ops.context_types_all),1);
num_trials_all = zeros(numel(app.ops.context_types_all),1);

hc_params.plot_dist_mat = 0;
hc_params.plot_clusters = 0;


for n_tt = 1:numel(app.ops.context_types_all)
    tt = app.ops.context_types_all(n_tt);
    tr_idx = tt == trial_types_wctx;
    
    if select_resp_cells
        cell_idx = stats1.cell_is_resp(:,n_tt);
    else
        cell_idx = true(stats1.num_cells,1);
    end
    
    num_tr = sum(tr_idx);
    num_cells = sum(cell_idx);
    
    data2 = trial_data_sort_wctx(cell_idx,:,tr_idx);
    data2_2d = reshape(data2,num_cells,[]);
    
    data_dim_est = f_ensemble_comp_data_dim2(data2_2d, est_params_pca);
    
    dim_est_tr(n_tt) = data_dim_est.dimensionality_corr;
    num_cells_all(n_tt) = num_cells;
    num_trials_all(n_tt) = num_tr;
    
    hc_out_cell = f_hcluster_wrap(data2_2d, hc_params);

    data2_2d_tr = reshape(data2,[], num_tr);
    hc_out_tr = f_hcluster_wrap(data2_2d_tr', hc_params);
    
    data2_2d_tr_sort = reshape(data2(hc_out_cell.dend_order,:,hc_out_tr.dend_order),num_cells,[]);
    dist1 = f_pdist_YS(data2_2d_tr', 'cosine');
    
    figure; 
    subplot(2,3,1);
    imagesc(data2_2d(hc_out_cell.dend_order,:));
    title(sprintf('Trial %d', tt));
    subplot(2,3,2);
    imagesc(1 - hc_out_cell.dist)
    sp = gca;
    sp.YDir = 'reverse';
    subplot(2,3,3);
    imagesc(data2_2d_tr_sort);
    subplot(2,3,4);
    imagesc(1 - dist1)
    sp = gca;
    sp.YDir = 'reverse';
    subplot(2,3,6);
    imagesc(1 - hc_out_tr.dist)
    sp = gca;
    sp.YDir = 'reverse';
    
end

figure; 
subplot(2,1,1);
bar(dim_est_tr)
subplot(2,1,2);
bar(num_cells_all)

end