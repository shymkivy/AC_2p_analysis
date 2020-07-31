function f_hclust_raster(trial_data_sort, trial_peaks, trial_types, sp, params, ops)

method = ops.dred_params.hclust.method;

[num_cells, num_bins, num_trials] = size(trial_data_sort);

if ops.dred_params.hclust.sort_raster
    dend_order1 = params.hclust_out_tr.dend_order;
else
    dend_order1 = 1:num_trials;
end
trial_data_sort2 = trial_data_sort(:,:,dend_order1);

trial_data_sort2 = reshape(trial_data_sort2, num_cells, []);

dend_order2 = params.hclust_out_cell.dend_order;

%[dend_order2, ~] = f_hcluster(trial_peaks, method);

subplot(sp)
imagesc(reshape(trial_data_sort2(dend_order2,:,:), num_cells, []));
title(sprintf('d%d', params.n_dset));

f_plot_trial_indicator(trial_types, dend_order1, num_bins, num_cells, ops)

end