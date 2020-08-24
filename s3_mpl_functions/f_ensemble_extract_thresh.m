function ens_out = f_ensemble_extract_thresh(coeffs, scores, num_clust, params, ops)
ensamble_method = f_get_param(params, 'ensamble_method', 'nmf');
cluster_method = f_get_param(params, 'cluster_method', 'hclust');    % 'hclust' or 'gmm'
cluster_method_cell = f_get_param(params, 'cluster_method_cell', 'hclust');
plot_stuff = f_get_param(params, 'plot_stuff', 0);

num_comps = size(scores,1);

%% first detect cells

% 
% figure;
% for n_cm = 1:(num_clust-1)
%     subplot((num_clust-1),1,n_cm);
%     stem(scores(n_cm,:))
% end


for n_comp = 1:(num_clust-1)
    subplot((num_clust-1),1,n_comp);
    
end

%% for each comp need thresh
X = coeffs';
comps_cell = if_get_comp_thresh(X);
 
X = scores;
comps_tr = if_get_comp_thresh(X);

figure;
stem(coeffs(:,1))

figure;
stem(scores(1,:))

figure; imagesc(coeffs*scores)

end

function comps_out = if_get_comp_thresh(X)

comps_out = cell((num_clust-1),1);

figure;
for n_comp = 1:(num_clust-1)
    factors = X(n_comp,:);
    center1 = median(factors);
    z_fac = sqrt(sum((factors-center1).^2)/(numel(factors)-1));
    %rms(coeffs(:,n_comp)-center1)
    
    zFactors = factors/z_fac;
    
    subplot((num_clust-1),1,n_comp); hold on;
    stem(zFactors)
    plot(ones(numel(factors),1)*2, '--r');
    
    comps_out.zFactors = zFactors;
    comps_out.z_thresh
    comps_out.sig_comps
    comps_out.insignif_comps
    %figure; histogram(coeffs(:,n_comp))
end

end
