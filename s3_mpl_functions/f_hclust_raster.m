function f_hclust_raster(trial_data_sort, trial_peaks, sp, params, ops)

method = ops.dred_params.hclust.method;

[num_cells, ~, ~] = size(trial_data_sort);


[dend_order, ~] = f_hcluster(trial_peaks, method);

subplot(sp)
imagesc(reshape(trial_data_sort(dend_order,:,:), num_cells, []));
title(sprintf('d%d', params.n_dset));



end