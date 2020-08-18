function ens_out = f_ensemble_extract_thresh(coeffs, scores, num_clust, params, ops)
ensamble_method = f_get_param(params, 'ensamble_method', 'nmf');
cluster_method = f_get_param(params, 'cluster_method', 'hclust');    % 'hclust' or 'gmm'
cluster_method_cell = f_get_param(params, 'cluster_method_cell', 'hclust');
plot_stuff = f_get_param(params, 'plot_stuff', 0);

num_comps = size(scores,1);

%% first detect cells

figure;
stem(coeffs(:,1))

figure;
stem(scores(1,:))

figure; imagesc(coeffs*scores)

end