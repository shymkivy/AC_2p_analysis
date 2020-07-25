function f_hclust_raster(trial_data_sort, sp, params)

[num_cells, num_bins, num_tr] = size(trial_data_sort);

subplot(sp)
imagesc(reshape(trial_data_sort, num_cells, []));
title(sprintf('d%d', params.n_dset));








end