function hclust_out = f_hcluster_trial3(X, params)
if ~exist('params', 'var')
    params = struct();
end

method = params.method;
metric = params.metric;
num_cells = size(X,1);

%% estimate num clust
if ~isfield(params, 'num_clust')
    if num_cells > 5
        k_list = 1:6;
    else
        k_list = 1:size(X,1);
    end

    E = evalclusters(X,'linkage','silhouette','klist',k_list, 'Distance', metric);
    num_clust = E.OptimalK;
else
    num_clust = params.num_clust;
end


%%

% warning is because some trial bins are nearly zero
[dend_order, clust_ident] = f_hcluster(X, method, num_clust);

if params.plot_sm
    image_Z = 1-squareform(pdist(X(dend_order,:), metric));
    figure;hold on; %subplot(sp); 
    imagesc(image_Z);
    %axis image;
    title(sprintf('h clust %s %s', method, metric));
    caxis([0 1]);
    axis tight;
    axis equal;
    xlabel('Trials');
    ylabel('Trials');
    %colorbar;
    %clim1 = caxis;
    sp = gca;
    sp.YDir = 'reverse';
end

%% add trial indicator

%f_plot_trial_indicator(trial_types, dend_order, 1, numel(trial_types), ops);

%imagesc(num_trials+(1:col_width),1:num_trials,permute(repmat(color_seq_tt,col_width,1,1),[2,1,3]));
%imagesc(num_trials+col_width+(1:col_width),1:num_trials,permute(repmat(color_seq_temporal,col_width,1,1),[2,1,3]));

%% plot clusters
if params.plot_sm
    if num_clust > 1
        ord1 = clust_ident(dend_order);
        for n_clust = 1:num_clust
            temp_list = find(ord1 == n_clust);
            subplot(sp); hold on;
            rectangle('Position',[temp_list(1)-0.5 temp_list(1)-0.5 numel(temp_list)-1+1 numel(temp_list)-1+1], 'EdgeColor', 'r','LineWidth',2);
        end
    end
end

hclust_out.num_clust = num_clust;
hclust_out.dend_order = dend_order;
hclust_out.clust_ident = clust_ident;
%hclust_out.clim = clim1;

end
