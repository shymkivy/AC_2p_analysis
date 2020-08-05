function f_mpl_trial_trial_dist(data, ops)

tn_to_dred = zeros(numel(ops.dred_params.trial_types_for_dist),1);
trial_type_tag = cell(numel(ops.dred_params.trial_types_for_dist),1);

for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data.(cond_name);
    
    dr_params.cond_name = cond_name;
    dr_params.colors_clust = cat(2,ops.colors_list,ops.colors_list,ops.colors_list,ops.colors_list,ops.colors_list);

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
        
    end
end





end