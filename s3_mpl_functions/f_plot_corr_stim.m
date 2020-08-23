function f_plot_corr_stim(trial_ave_vec, tt_ind, stim_ind, ops)

corr_list = cell(numel(ops.regions_to_analyze),1);
for n_cond = 1:numel(ops.regions_to_analyze)
    num_dsets = size(trial_ave_vec{n_cond},2);
    corr_list{n_cond} = zeros(num_dsets,1);
    corr_list2 = zeros(size(stim_ind,1), num_dsets);
    for n_stim = 1:size(stim_ind,1)
        
        for n_dset = 1:num_dsets
            temp_data = trial_ave_vec{n_cond}{tt_ind,n_dset};
            corr_list2(n_stim,n_dset) = 1-pdist2(temp_data(:,stim_ind(n_stim,1))',temp_data(:,stim_ind(n_stim,2))','cosine');
        end
        
    end
    corr_list{n_cond} = corr_list2(:);
end

p_val_mat = f_get_tt_stats(corr_list);
figure; imagesc(p_val_mat);
f_plot_dset_deets(corr_list, ops)


end