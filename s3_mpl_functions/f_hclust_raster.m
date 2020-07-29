function f_hclust_raster(trial_data_sort, trial_peaks, sp, params, ops)

method = ops.dred_params.hclust.method;

[num_cells, ~, ~] = size(trial_data_sort);

if ops.dred_params.hclust.sort_raster
    trial_data_sort2 = trial_data_sort(:,:,params.hclust_out.dend_order);
end

trial_data_sort2 = reshape(trial_data_sort2, num_cells, []);

[dend_order, ~] = f_hcluster(trial_peaks, method);



subplot(sp)
imagesc(reshape(trial_data_sort2(dend_order,:,:), num_cells, []));
title(sprintf('d%d', params.n_dset));



end