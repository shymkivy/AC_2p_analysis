function f_mpl_trial_trial_dist(data, ops)

tn_to_dred = zeros(numel(ops.dred_params.trial_types_for_dist),1);
trial_type_tag = cell(numel(ops.dred_params.trial_types_for_dist),1);

for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data.(cond_name);
    
    dr_params.cond_name = cond_name;
    dr_params.colors_clust = cat(2,ops.colors_list,ops.colors_list,ops.colors_list,ops.colors_list,ops.colors_list);
    dr_params.num_clust = ops.dred_params.hclust.num_clust{1};
    
    
    for n_dset = 1:cdata.num_dsets
        disp([cond_name, ' dset ' num2str(n_dset)]);

        trial_types = cdata.trial_types_pr{n_dset};
        trial_peaks = cdata.tuning_all{n_dset}.peak_tuning_full_resp.fr_peak_mag;
        trial_data_sort_sm_pr = cdata.trial_data_sort_sm_pr{n_dset};
        
        [tn_to_dred(1), trial_type_tag{1}] = f_select_trial_type(ops.dred_params.trial_types_for_dist{1}, cdata, n_dset, ops);
        [tn_to_dred(2), trial_type_tag{2}] = f_select_trial_type(ops.dred_params.trial_types_for_dist{2}, cdata, n_dset, ops);
        tt_to_dred = ops.context_types_all(tn_to_dred);
        
        trials_idx_dred = logical(sum(trial_types == tt_to_dred(:)' ,2));
        trial_peaks_dred = trial_peaks(:,trials_idx_dred);
        trial_types_dred = trial_types(trials_idx_dred);
        trial_data_sort_sm_pr = trial_data_sort_sm_pr(:,:,trials_idx_dred);
        
        if ops.dred_params.use_responsive_cells
            resp_cells = cdata.peak_tuned_trials_combined{n_dset};
            resp_cells = logical(sum(resp_cells(:,tn_to_dred),2));
            trial_peaks_dred = trial_peaks_dred(resp_cells,:);
            trial_data_sort_sm_pr = trial_data_sort_sm_pr(resp_cells,:,:);
        end
        
        params2.method = 'ward'; % cosine, ward
        params2.metric = 'euclidean'; % cosine squaredeuclidean
        params2.plot_sm = 1;
        hclust_out = f_hcluster_trial3(trial_peaks_dred', params2);
        
        %sp_h_tr = figure;
        %sp_h_tr{n_dset} = subplot(3,5,n_dset);
        hclust_out_tr = f_hcluster_trial2(trial_peaks_dred, trial_types_dred, [], dr_params, ops);
        %dr_params.hclust_out_tr = hclust_out_tr{n_dset};
        
        %% hclustering cells 
        figure;
        sp_h_cell{n_dset} = subplot(3,5,n_dset);
        hclust_out_cell{n_dset} = f_hcluster_cell(trial_peaks_dred, trial_types_dred, sp_h_cell{n_dset}, dr_params, ops);
        dr_params.hclust_out_cell = hclust_out_cell{n_dset};
        %%
        %f_tsne(trial_peaks)

        %%
        %ops.dred_params.hclust.sort_raster = 1;
        figure;
        sp_ras{n_dset} = subplot(3,5,n_dset);
        f_hclust_raster(trial_data_sort_sm_pr, trial_peaks_dred, trial_types_dred, sp_ras{n_dset}, dr_params, ops);
        
        
    end
end





end