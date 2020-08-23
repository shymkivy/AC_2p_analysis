function clust_out = f_get_clust_params(X, clust_in)

num_dred_comps = size(X,2);
clust_label = clust_in.clust_label;
num_clust = numel(clust_in.clust_label);

clust_centers = zeros(num_clust,num_dred_comps);
clust_mag = zeros(num_clust,1);
num_cells = zeros(num_clust,1);
for n_clust = 1:num_clust
    clust1 = clust_label(n_clust);
    clust_centers(n_clust,:) = mean(X(clust_in.clust_ident == clust1,:));
    clust_mag(n_clust) = norm(clust_centers(n_clust,:));
    num_cells(n_clust) = sum(clust_in.clust_ident == clust1);
end


clust_out = clust_in;
clust_out.clust_centers = clust_centers;
clust_out.clust_mag = clust_mag;
clust_out.num_cells = num_cells;
end