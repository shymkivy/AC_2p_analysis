function [dend_order, clust_ident] = f_hierarch_clust(data, dend_thresh)

if ~exist('dend_thresh', 'var')
    dend_thresh = 50;
end
method1 = 'ward';

Z = linkage(data,method1);
%Z = linkage(1-squareform(pdist(data,'cosine')),method1);

f1 = figure;
[~, ~, dend_order] = dendrogram(Z, 1000,'ColorThreshold',dend_thresh,'Orientation','left');
close(f1);

f2 = if_plot_trial_trial_image(1-squareform(pdist(data(dend_order,:), 'cosine')), 'cosine');
if_plot_trial_trial_image(1-squareform(pdist(data(dend_order,:), 'euclidean')), 'euclidean');
%if_plot_trial_trial_image(squareform(pdist(data(dend_order,:), 'minkowski',1)), 'minkowski 1');
%if_plot_trial_trial_image(squareform(pdist(data(dend_order,:), 'minkowski',2)), 'minkowski 2');


num_clust = input('How many clusters?');

clust_ident = cluster(Z, 'MaxClust', num_clust);

% if num_clust > 1
%     figure(f1);
%     [~, clust_ident, ~] = dendrogram(linkage(data,method1), num_clust,'ColorThreshold',1000,'Orientation','left');
% else
%     clust_ident = ones(size(data,1),1);
% end

ord1 = clust_ident(dend_order);
figure(f2); hold on;
for n_clust = 1:num_clust
    temp_list = find(ord1 == n_clust);
    rectangle('Position',[temp_list(1) temp_list(1) numel(temp_list)-1 numel(temp_list)-1], 'EdgeColor', 'r','LineWidth',2)
end

end

function fhandle = if_plot_trial_trial_image(image_Z, metric)


fhandle = figure;
imagesc(image_Z);
axis image;
title(sprintf('tr-tr similarity ward sorted, %s', metric));
%caxis([-.2 .8]);
axis equal tight;
xlabel('Trials');
ylabel('Trials');
colorbar;

end