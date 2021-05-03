function f_hclust_raster(trial_data_sort, trial_types, sp, params, ops)

ndims1 = ndims(trial_data_sort);
if ndims1 == 2
    [num_cells, num_trials] = size(trial_data_sort);
    num_bins = 1;
elseif ndims1 == 3
    [num_cells, num_bins, num_trials] = size(trial_data_sort);
end

if ops.dred_params.hclust.sort_raster
    dend_order_trials = params.dend_order_trials;
else
    dend_order_trials = 1:num_trials;
end

if ndims1 == 2
    trial_data_sort3 = trial_data_sort(:,dend_order_trials);
elseif ndims1 == 3
    trial_data_sort2 = trial_data_sort(:,:,dend_order_trials);
    trial_data_sort3 = reshape(trial_data_sort2, num_cells, []);
end

dend_order_cells = params.dend_order_cells;

%[dend_order2, ~] = f_hcluster(trial_peaks, method);

raster2 = reshape(trial_data_sort3(dend_order_cells,:,:), num_cells, []);

subplot(sp)
imagesc(raster2);
title(sprintf('d%d, corr=%.2f', params.n_dset, params.dim_corr));

trial_numbers = f_tt_to_tn(trial_types, ops);
added_width = f_plot_trial_indicator(trial_numbers, dend_order_trials, num_bins, num_cells, ops.context_types_all_colors2);

raster3 = zeros(size(raster2)+[added_width*2 0]);
clust_ident_trials = reshape(repmat(params.clust_ident_trials(dend_order_trials)',num_bins,1),[],1);

f_plot_cell_indicator(raster3, params.clust_ident_cells(dend_order_cells), ops);
f_plot_trial_indicator2(raster3, clust_ident_trials, 1, ops);

end