function clust_out = f_ensemble_clust_cell(coeffs, scores, num_ens, raster_norm, params)
ensamble_method = f_get_param(params, 'ensamble_method', 'nmf');
cluster_method = f_get_param(params, 'cluster_method', 'hclust');    % 'hclust' or 'gmm'
cluster_method_cell = f_get_param(params, 'cluster_method_cell', 'hclust');
plot_stuff = f_get_param(params, 'plot_stuff', 0);

num_comps = num_ens;
num_clust = num_ens + 1;

%%

dist_out = f_rbf_kernel(X, X);

[dend_order, clust_ident, Z] = f_hcluster(X, 'rbf', num_clust);

figure; imagesc(dist_out(dend_order,dend_order))
%%

X = coeffs;
%X = firing_rate_LR;



%% cluster cells
if strcmpi(cluster_method_cell, 'hclust')
    %% cluster with hclust
    params2.method = 'rbf'; % cosine, ward, rbf
    params2.metric = 'cosine'; % cosine squaredeuclidean
    params2.plot_dist_mat = plot_stuff;
    params2.plot_clusters = plot_stuff;
    params2.num_clust = num_clust;
    params2.XY_label = 'Cells';
    clust_out_cell = f_hcluster_wrap(X, params2);
    %gscatter(X(:,1),X(:,2),hclust_out.clust_ident);
    
elseif strcmpi(cluster_method_cell, 'gmm')
    %% clust with gmm, can find best regularizer
    optimize_reg_val = 0;
    if optimize_reg_val
        num_vals = 50;
        num_reps = 10;
        mean_acc = zeros(num_vals,num_reps);
        %rg_list = linspace(0.0001, 1, num_vals);
        rg_list = logspace(-4, -2, num_vals);
        for n_rg = 1:numel(rg_list)
            for n_rep = 1:num_reps
                params3.metric = 'cosine'; % cosine squaredeuclidean
                params3.RegularizationValue = rg_list(n_rg);
                params3.num_clust = num_comps+1;
                gmmclust_out = f_gmmcluster_trial(X, params3);
                %gscatter(X(:,1),X(:,2),gmmclust_out.clust_ident);
                eval_gmm = f_evaluate_ens_result(gmmclust_out.clust_ident, ens_list_gt, 0);
                mean_acc(n_rg, n_rep) = mean(eval_gmm.accuracy);
            end
        end
        f_rg = figure;
        shadedErrorBar(log10(rg_list), mean(mean_acc,2), std(mean_acc, [], 2)/sqrt(num_reps-1));
        xlabel('log10(rg val)')
        title('Click on optimal regression value');
        [x,~] = ginput(1);
        close(f_rg);
        params3.RegularizationValue = 10^x;
    else
        params3.RegularizationValue = 0.006;
    end
    params3.metric = 'sqEeuclidean'; % cosine squaredeuclidean
    params3.num_clust = num_clust;
    clust_out_cell = f_gmmcluster_trial(X, params3);
    [~, clust_out_cell.dend_order] = sort(clust_out_cell.clust_ident);
end

clust_params_cell = if_get_clust_params(X, clust_out_cell);

%proj_out = f_ens_project_clust(coeffs, scores, clust_params_cell, clust_params_tr, raster_norm);


end

function clust_params = if_get_clust_params(X, clust_out)

num_dred_comps = size(X,2);

clust_label = unique(clust_out.clust_ident);
num_clust = numel(clust_label);

clust_centers = zeros(num_clust,num_dred_comps);
clust_mag = zeros(num_clust,1);
cell_num = zeros(num_clust,1);
for n_clust = 1:num_clust
    clust1 = clust_label(n_clust);
    clust_centers(n_clust,:) = mean(X(clust_out.clust_ident == clust1,:));
    clust_mag(n_clust) = norm(clust_centers(n_clust,:));
    cell_num(n_clust) = sum(clust_out.clust_ident == clust1);
end

[~, ens_order] = sort(clust_mag);
[~, ens_order2] = sort(ens_order);
clust_ident2 = ens_order2(clust_out.clust_ident)-1;

ens_list = cell(num_clust-1,1);
for n_ens = 1:(num_clust-1)
    ens_list{n_ens} = find(clust_ident2 == (n_ens));
end

clust_params.clust_label = clust_label - 1;
clust_params.ens_list = ens_list;
clust_params.residual_list = find(clust_ident2 == 0);
clust_params.clust_ident = clust_ident2(:);
clust_params.dend_order = clust_out.dend_order(:);
end