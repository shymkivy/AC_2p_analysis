function f_hcluster_trial(trial_peaks, trial_types, method, num_clust, params)

if isempty(num_clust)
    num_clust = 1;
end


[dend_order, clust_ident] = f_hcluster(trial_peaks', method, num_clust);

figure;
sp{1} = subplot(1,3,1);
if_plot_trial_trial_image(1-squareform(pdist(trial_peaks(:,dend_order)', 'cosine')), 'cosine');
sp{2} = subplot(1,3,2);
if_plot_trial_trial_image(1-squareform(pdist(trial_peaks(:,dend_order)', 'euclidean')), 'euclidean');
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
suptitle(sprintf('%s dset %d; %s clust=%d; trials:[%s], %d cluster', params.cond_name, params.n_dset, method, num_clust, num2str(params.tt_to_dred(:)')));

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