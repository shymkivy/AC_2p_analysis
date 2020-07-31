function hclust_out = f_hcluster_cell(trial_peaks, trial_types, sp, params, ops)

num_clust = params.num_clust;
method = ops.dred_params.hclust.method;
metric = ops.dred_params.hclust.plot_metric;
num_cells = size(trial_peaks,1);

if isempty(num_clust)
    num_clust = 1;
end

% warning is because some trial bins are nearly zero
[dend_order, clust_ident] = f_hcluster(trial_peaks, method, num_clust);

%figure; imagesc(color_seq_temporal)

image_Z = 1-squareform(pdist(trial_peaks(dend_order,:), metric));
subplot(sp); hold on;
imagesc(image_Z);
%axis image;
title(sprintf('d%d, %d cells', params.n_dset,num_cells));
caxis([0 1]);
axis tight;
axis equal;
xlabel('Cells');
ylabel('Cells');
%colorbar;
clim1 = caxis;
sp.YDir = 'reverse';

%% add trial indicator

%f_plot_trial_indicator(trial_types, dend_order, 1, numel(trial_types), ops);

%imagesc(num_trials+(1:col_width),1:num_trials,permute(repmat(color_seq_tt,col_width,1,1),[2,1,3]));
%imagesc(num_trials+col_width+(1:col_width),1:num_trials,permute(repmat(color_seq_temporal,col_width,1,1),[2,1,3]));

%% plot clusters
if num_clust > 1
    ord1 = clust_ident(dend_order);
    for n_clust = 1:num_clust
        temp_list = find(ord1 == n_clust);
        subplot(sp); hold on;
        rectangle('Position',[temp_list(1)-0.5 temp_list(1)-0.5 numel(temp_list)-1+1 numel(temp_list)-1+1], 'EdgeColor', 'r','LineWidth',2);
    end
end

hclust_out.dend_order = dend_order;
hclust_out.clust_ident = clust_ident;
hclust_out.clim = clim1;

end
