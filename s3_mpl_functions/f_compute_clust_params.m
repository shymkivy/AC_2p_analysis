function clust_params = f_compute_clust_params(pk_mag, clust_ident)

clust_name = unique(clust_ident);
num_clust = numel(clust_name);            
            
clust_cent = zeros(num_clust,3);
clust_num_cells = zeros(num_clust,1);
clust_mean_dist = zeros(num_clust,1);
for n_cl = 1:num_clust
    clust_cells = clust_ident == clust_name(n_cl);
    clust_num_cells(n_cl) = sum(clust_cells);
    clust_cent(n_cl,:) = mean(pk_mag(clust_cells,:));
    clust_mean_dist(n_cl) = mean(sqrt(sum((pk_mag(clust_cells,:) - clust_cent(n_cl,:)).^2,2)));
end

clust_params.clust_name = clust_name;
clust_params.clust_cent = clust_cent;
clust_params.clust_num_cells = clust_num_cells;
clust_params.clust_mean_dist = clust_mean_dist;
end