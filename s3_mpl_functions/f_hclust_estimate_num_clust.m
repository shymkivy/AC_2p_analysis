function est = f_hclust_estimate_num_clust(trial_peaks)

%method = ops.dred_params.hclust.method;
plot_metric = 'cosine';
make_plot = 1;

if size(trial_peaks,1) > 6
    k_list = 1:6;
else
    k_list = 1:size(trial_peaks,1);
end



E1 = evalclusters(trial_peaks','linkage','silhouette','klist',k_list, 'Distance', plot_metric);
E2 = evalclusters(trial_peaks','linkage','CalinskiHarabasz','klist',k_list);
E3 = evalclusters(trial_peaks','linkage','DaviesBouldin','klist',k_list);
E4 = evalclusters(trial_peaks','linkage','gap','klist',k_list, 'Distance', plot_metric);

if make_plot
    figure;
    subplot(2,2,1); plot(E1);
    title(sprintf('OptimalK = %d', E1.OptimalK))
    subplot(2,2,2); plot(E2);
    title(sprintf('OptimalK = %d', E2.OptimalK))
    subplot(2,2,3); plot(E3);
    title(sprintf('OptimalK = %d', E3.OptimalK))
    subplot(2,2,4); plot(E4);
    title(sprintf('OptimalK = %d', E4.OptimalK))
    %suptitle(sprintf('%s dset %d; %s most likely clust num; trials:[%s]', params.cond_name, params.n_dset, method, num2str(params.tt_to_dred(:)')));
end

est.silhouette = E1;
est.CalinskiHarabasz = E2;
est.DaviesBouldin = E3;
est.gap = E4;

end
