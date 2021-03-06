function f_hcluster_trial(trial_peaks, trial_types, method, num_clust, params, ops)

if isempty(num_clust)
    num_clust = 1;
end

if size(trial_peaks,1) > 9
    k_list = 1:10;
else
    k_list = 1:size(trial_peaks,1);
end

figure; 
E = evalclusters(trial_peaks','linkage','silhouette','klist',k_list, 'Distance', 'cosine');
subplot(2,2,1); plot(E);
E = evalclusters(trial_peaks','linkage','CalinskiHarabasz','klist',k_list);
subplot(2,2,2); plot(E);
E = evalclusters(trial_peaks','linkage','DaviesBouldin','klist',k_list);
subplot(2,2,3); plot(E);
E = evalclusters(trial_peaks','linkage','gap','klist',k_list, 'Distance', 'cosine');
subplot(2,2,4); plot(E);
suptitle(sprintf('%s dset %d; %s most likely clust num; trials:[%s]', params.cond_name, params.n_dset, method, num2str(params.tt_to_dred(:)')));


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

figure;
sp{1} = subplot(1,3,1); hold on;
if_plot_trial_trial_image(1-squareform(pdist(trial_peaks(:,dend_order)', 'cosine')), 'cosine');
sp{1}.YDir = 'reverse';
imagesc(1:num_trials,num_trials+(1:col_width),repmat(color_seq_tt,col_width,1,1));
imagesc(1:num_trials,num_trials+col_width+(1:col_width),repmat(color_seq_temporal,col_width,1,1));
%imagesc(num_trials+(1:col_width),1:num_trials,permute(repmat(color_seq_tt,col_width,1,1),[2,1,3]));
%imagesc(num_trials+col_width+(1:col_width),1:num_trials,permute(repmat(color_seq_temporal,col_width,1,1),[2,1,3]));


sp{2} = subplot(1,3,2); hold on;
if_plot_trial_trial_image(1-squareform(pdist(trial_peaks(:,dend_order)', 'euclidean')), 'euclidean');
sp{2}.YDir = 'reverse';
imagesc(1:num_trials,num_trials+(1:col_width),repmat(color_seq_tt,col_width,1,1));
imagesc(1:num_trials,num_trials+col_width+(1:col_width),repmat(color_seq_temporal,col_width,1,1));
%imagesc(num_trials+(1:col_width),1:num_trials,permute(repmat(color_seq_tt,col_width,1,1),[2,1,3]));
%imagesc(num_trials+col_width+(1:col_width),1:num_trials,permute(repmat(color_seq_temporal,col_width,1,1),[2,1,3]));

%if_plot_trial_trial_image(squareform(pdist(data(dend_order,:), 'minkowski',1)), 'minkowski 1');
%if_plot_trial_trial_image(squareform(pdist(data(dend_order,:), 'minkowski',2)), 'minkowski 2');



% plot clusters
if num_clust > 1
    ord1 = clust_ident(dend_order);
    for n_clust = 1:num_clust
        temp_list = find(ord1 == n_clust);
        for n_sp = 1:numel(sp)
            subplot(sp{n_sp}); hold on;
            rectangle('Position',[temp_list(1)-0.5 temp_list(1)-0.5 numel(temp_list)-1+1 numel(temp_list)-1+1], 'EdgeColor', 'r','LineWidth',2);
        end
    end
end

Y = tsne(trial_peaks');
subplot(1,3,3); hold on;
%gscatter(Y(:,1),Y(:,2),clust_ident);
leg1 = cell(num_clust,1);
for n_cl = 1:num_clust
    scatter(Y(clust_ident==n_cl,1),Y(clust_ident==n_cl,2),25,params.colors_clust{n_cl}, 'filled');
    leg1{n_cl} = num2str(n_cl);
end
axis tight;
axis square;
title('T-SNE');
if num_clust > 1
    legend(leg1);
end
suptitle(sprintf('%s dset %d; %s clust=%d; trials:[%s]', params.cond_name, params.n_dset, method, num_clust, num2str(params.tt_to_dred(:)')));

end

function clim1 = if_plot_trial_trial_image(image_Z, metric)
imagesc(image_Z);
%axis image;
title(sprintf('tr-tr sim sort, %s', metric));
%caxis([-.2 .8]);
axis tight;
axis equal;
xlabel('Trials');
ylabel('Trials');
%colorbar;
clim1 = caxis;
end