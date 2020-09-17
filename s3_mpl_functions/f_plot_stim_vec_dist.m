function f_plot_stim_vec_dist(trial_ave_vec, tt_ind, stim_ind, ops, dist_metric)

corr_list = cell(numel(ops.regions_to_analyze),1);
for n_cond = 1:numel(ops.regions_to_analyze)
    num_dsets = size(trial_ave_vec{n_cond},2);
    corr_list{n_cond} = zeros(num_dsets,1);
    corr_list2 = zeros(size(stim_ind,1), num_dsets);
    for n_stim = 1:size(stim_ind,1)
        
        for n_dset = 1:num_dsets
            temp_data = trial_ave_vec{n_cond}{tt_ind,n_dset};
            
            if strcmpi(dist_metric, 'cosineSI')
                corr_list2(n_stim,n_dset) = 1-pdist2(temp_data(:,stim_ind(n_stim,1))',temp_data(:,stim_ind(n_stim,2))',dist_metric);
            else
                corr_list2(n_stim,n_dset) = pdist2(temp_data(:,stim_ind(n_stim,1))',temp_data(:,stim_ind(n_stim,2))',dist_metric);
            end
        end
        
    end
    corr_list{n_cond} = corr_list2(:);
end

p_val_mat = f_get_tt_stats(corr_list);
figure;
subplot(2,1,2);
imagesc(p_val_mat);
caxis([0 1]);
sp1 = subplot(2,1,1);
f_plot_dset_deets(corr_list, ops, sp1)





end