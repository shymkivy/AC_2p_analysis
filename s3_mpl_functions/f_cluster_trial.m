function f_cluster_trial(trial_peaks, trial_types, params)

%trial_types_cat = unique(trial_types);

figure;
sp{1} = subplot(1,3,1);
sp{2} = subplot(1,3,2);
[~, clust_ident] = f_hierarch_clust(trial_peaks', 2, sp);

Y = tsne(trial_peaks');
subplot(1,3,3);
gscatter(Y(:,1),Y(:,2),clust_ident);
axis tight;
axis equal;
title('T-SNE');
suptitle(sprintf('%s dset %d, trials:[%s]', params.cond_name, params.n_dset,num2str(params.tn_to_dred)));

end