function f_hcluster_cond(cdata, dr_params, ops)


hclust_out = cell(cdata.num_dsets,1);
fig_h = figure;
sp_h = cell(cdata.num_dsets,1);
fig_ras = figure;
sp_ras = cell(cdata.num_dsets,1);

for n_dset = 1:cdata.num_dsets
    disp([dr_params.cond_name, ' dset ' num2str(n_dset)]);

    trial_types = cdata.trial_types_pr{n_dset};
    trial_peaks = cdata.tuning_all{n_dset}.peak_tuning_full_resp.fr_peak_mag;
    trial_data_sort_sm_pr = cdata.trial_data_sort_sm_pr{n_dset};
    
    %% select trials
    [tn_to_dred, trial_type_tag] = f_select_trial_type(dr_params.tt_to_dred_input, cdata, n_dset, ops);
    tt_to_dred = ops.context_types_all(tn_to_dred);
    
    trials_idx_dred = logical(sum(trial_types == tt_to_dred(:)' ,2));
    trial_peaks_dred = trial_peaks(:,trials_idx_dred);
    trial_types_dred = trial_types(trials_idx_dred);
    trial_data_sort_sm_pr = trial_data_sort_sm_pr(:,:,trials_idx_dred);
    
    %% dset specific params
    dr_params.n_dset = n_dset;
    dr_params.trial_type_tag = trial_type_tag;
    dr_params.tn_to_dred = tn_to_dred;
    dr_params.tt_to_dred = tt_to_dred;
    dr_params.volume_period = cdata.proc_data{n_dset}.frame_data.volume_period;
    dr_params.trial_t = cdata.trial_window_t{n_dset};
    dr_params.ctx_mmn = ops.context_types_all(cdata.ctx_mmn{n_dset});
    

    %% select responsive cells
    if ops.dred_params.use_responsive_cells
        resp_cells = cdata.peak_tuned_trials_combined{n_dset};
        resp_cells = logical(sum(resp_cells(:,tn_to_dred),2));
        trial_peaks_dred = trial_peaks_dred(resp_cells,:);
        trial_data_sort_sm_pr = trial_data_sort_sm_pr(resp_cells,:,:);
    end
    
    if sum(resp_cells)>1
        %%
        %f_hclust_estimate_num_clust(trial_peaks_dred, dr_params, ops)

        %% hclustering
        figure(fig_h);
        sp_h{n_dset} = subplot(3,5,n_dset);
        hclust_out{n_dset} = f_hcluster_trial2(trial_peaks_dred, trial_types_dred, sp_h{n_dset}, dr_params, ops);
        dr_params.hclust_out = hclust_out{n_dset};
        %%
        %f_tsne(trial_peaks)

        %%
        figure(fig_ras);
        sp_ras{n_dset} = subplot(3,5,n_dset);
        f_hclust_raster(trial_data_sort_sm_pr, trial_peaks_dred, sp_ras{n_dset}, dr_params, ops);

    end
end
figure(fig_h);
suptitle(sprintf('%s; %s clust=%d; trials:[%s]', dr_params.cond_name, ops.dred_params.hclust.method, dr_params.num_clust, num2str(tt_to_dred(:)')));
figure(fig_ras);
suptitle(sprintf('%s; %s clust=%d; trials:[%s]', dr_params.cond_name, ops.dred_params.hclust.method, dr_params.num_clust, num2str(tt_to_dred(:)')));


end