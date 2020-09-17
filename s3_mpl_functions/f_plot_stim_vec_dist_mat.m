function f_plot_stim_vec_dist_mat(trial_ave_vec, tt_tag, ops)

corr_all = cell(numel(ops.regions_to_analyze), numel(ops.dred_params.trial_types_for_dist));
for n_tt = 1:numel(ops.dred_params.trial_types_for_dist)
    for n_cond = 1:numel(ops.regions_to_analyze)
        temp_data_full = cat(1,trial_ave_vec{n_cond}{n_tt,:});
        %corr_sim = corr(temp_data_full,temp_data_full);
        corr_sim = 1-pdist2(temp_data_full',temp_data_full','cosine');
        corr_all{n_cond, n_tt} = corr_sim;       
    end
    %mean_tt = mean(cat(3,corr_all{:,n_tt}),3);
    for n_cond = 1:numel(ops.regions_to_analyze)
        cond_name = ops.regions_to_analyze{n_cond};
        figure; imagesc(corr_all{n_cond, n_tt});
        %caxis([-.3 .3])
        title(sprintf('%s, %s, corr',tt_tag{n_tt},  cond_name));
    end
end


end