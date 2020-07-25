function f_hclust_estimate_num_clust(trial_peaks, params, ops)

method = ops.dred_params.hclust.method;

if size(trial_peaks,1) > 9
    k_list = 1:10;
else
    k_list = 1:size(trial_peaks,1);
end

figure; 
E = evalclusters(trial_peaks','linkage','silhouette','klist',k_list, 'Distance', params.plot_metric);
subplot(2,2,1); plot(E);
E = evalclusters(trial_peaks','linkage','CalinskiHarabasz','klist',k_list);
subplot(2,2,2); plot(E);
E = evalclusters(trial_peaks','linkage','DaviesBouldin','klist',k_list);
subplot(2,2,3); plot(E);
E = evalclusters(trial_peaks','linkage','gap','klist',k_list, 'Distance', params.plot_metric);
subplot(2,2,4); plot(E);
suptitle(sprintf('%s dset %d; %s most likely clust num; trials:[%s]', params.cond_name, params.n_dset, method, num2str(params.tt_to_dred(:)')));








end
