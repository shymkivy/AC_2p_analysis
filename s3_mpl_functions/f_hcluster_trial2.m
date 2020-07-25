function hclust_out = f_hcluster_trial2(trial_peaks, trial_types, sp, params, ops)

num_clust = params.num_clust;
method = ops.dred_params.hclust.method;
metric = ops.dred_params.hclust.plot_metric;
num_cells = size(trial_peaks,1);

if isempty(num_clust)
    num_clust = 1;
end

[dend_order, clust_ident] = f_hcluster(trial_peaks', method, num_clust);
num_trials = numel(trial_types);

trial_order = 1:num_trials;
trial_types_sort = trial_types(dend_order);
trial_order_sort = trial_order(dend_order);
gray_cmap = repmat(linspace(0.2,1,num_trials),3,1);

color_seq_tt = zeros(1,numel(trial_types),3);
color_seq_temporal = zeros(1,numel(trial_types),3);
for n_tr = 1:num_trials
    color_seq_tt(1,n_tr,:) = ops.context_types_all_colors(trial_types_sort(n_tr) == ops.context_types_all,:,:);
    color_seq_temporal(1,n_tr,:) = gray_cmap(:,trial_order_sort(n_tr));
end
col_width = ceil(num_trials/50);

%figure; imagesc(color_seq_temporal)

image_Z = 1-squareform(pdist(trial_peaks(:,dend_order)', metric));

subplot(sp); hold on;

imagesc(image_Z);
%axis image;
title(sprintf('d%d, %d cells', params.n_dset,num_cells));
caxis([0 1]);
axis tight;
axis equal;
xlabel('Trials');
ylabel('Trials');
%colorbar;
clim1 = caxis;

sp.YDir = 'reverse';
imagesc(1:num_trials,num_trials+(1:col_width),repmat(color_seq_tt,col_width,1,1));
imagesc(1:num_trials,num_trials+col_width+(1:col_width),repmat(color_seq_temporal,col_width,1,1));
%imagesc(num_trials+(1:col_width),1:num_trials,permute(repmat(color_seq_tt,col_width,1,1),[2,1,3]));
%imagesc(num_trials+col_width+(1:col_width),1:num_trials,permute(repmat(color_seq_temporal,col_width,1,1),[2,1,3]));

% plot clusters
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
