function hclust_out = f_hcluster_cell(trial_peaks, cell_types, params, ops)
n_dset = f_get_param(params, 'n_dset', 0);
num_clust = f_get_param(params, 'num_clust');
method = f_get_param(params, 'method', 'cosine');
metric = f_get_param(params, 'metric', 'cosine');
sp = f_get_param(params, 'subplot_ptr');
plot_dist_mat = f_get_param(params, 'plot_dist_mat', 1);
plot_clusters = f_get_param(params, 'plot_clusters', 1);

num_cells = size(trial_peaks,1);
if isempty(num_clust)
    num_clust = 1;
end

% warning is because some trial bins are nearly zero
[dend_order, clust_ident, Z] = f_hcluster(trial_peaks, method, num_clust);

%figure; imagesc(color_seq_temporal)

dist1 = pdist(trial_peaks(dend_order,:), metric);

hclust_out.dist = dist1;
hclust_out.dend_order = dend_order;
hclust_out.clust_ident = clust_ident;
hclust_out.Z = Z;

if plot_dist_mat
    image_Z = 1-squareform(dist1);
    if isempty(sp)
        figure;
        sp = gca;
    else
        subplot(sp);
    end
    hold on;
    imagesc(image_Z);
    %axis image;
    title(sprintf('d%d, %d cells', n_dset,num_cells));
    caxis([0 1]);
    axis tight;
    axis equal;
    xlabel('Cells');
    ylabel('Cells');
    %colorbar;
    clim1 = caxis;
    sp.YDir = 'reverse';

    %% add trial indicator
    if ~isempty(cell_types)
        f_plot_trial_indicator(cell_types, dend_order, 1, numel(cell_types), ops);
    end

    %imagesc(num_trials+(1:col_width),1:num_trials,permute(repmat(color_seq_tt,col_width,1,1),[2,1,3]));
    %imagesc(num_trials+col_width+(1:col_width),1:num_trials,permute(repmat(color_seq_temporal,col_width,1,1),[2,1,3]));

    %% plot clusters

    if plot_clusters
        ord1 = clust_ident(dend_order);
        for n_clust = 1:num_clust
            temp_list = find(ord1 == n_clust);
            subplot(sp); hold on;
            rectangle('Position',[temp_list(1)-0.5 temp_list(1)-0.5 numel(temp_list)-1+1 numel(temp_list)-1+1], 'EdgeColor', 'r','LineWidth',2);
        end
    end
    hclust_out.clim = clim1;
end

end
