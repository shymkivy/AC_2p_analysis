function [dend_order, clust_ident] = f_hierarch_clust(data, num_clust, fig_hand)

dend_thresh = 50;

if ~exist('fig_hand', 'var')
    fig_hand{1} = figure;
    fig_hand{2} = figure;
end

%Z = linkage(data,'ward');
Z = linkage(pdist(data,'cosine'),'average');

f1 = figure;
[~, ~, dend_order] = dendrogram(Z, 1000,'ColorThreshold',dend_thresh,'Orientation','left');
close(f1);

subplot(fig_hand{1});
if_plot_trial_trial_image(1-squareform(pdist(data(dend_order,:), 'cosine')), 'cosine');
subplot(fig_hand{2});
if_plot_trial_trial_image(1-squareform(pdist(data(dend_order,:), 'euclidean')), 'euclidean');
%if_plot_trial_trial_image(squareform(pdist(data(dend_order,:), 'minkowski',1)), 'minkowski 1');
%if_plot_trial_trial_image(squareform(pdist(data(dend_order,:), 'minkowski',2)), 'minkowski 2');

if ~exist('num_clust', 'var')
    num_clust = input('How many clusters?');
end


clust_ident = cluster(Z, 'MaxClust', num_clust);

% if num_clust > 1
%     figure(f1);
%     [~, clust_ident, ~] = dendrogram(linkage(data,method1), num_clust,'ColorThreshold',1000,'Orientation','left');
% else
%     clust_ident = ones(size(data,1),1);
% end

ord1 = clust_ident(dend_order);

for n_clust = 1:num_clust
    temp_list = find(ord1 == n_clust);
    subplot(fig_hand{1}); hold on;
    rectangle('Position',[temp_list(1)-0.5 temp_list(1)-0.5 numel(temp_list)-1+1 numel(temp_list)-1+1], 'EdgeColor', 'r','LineWidth',2);
    subplot(fig_hand{2}); hold on;
    rectangle('Position',[temp_list(1)-0.5 temp_list(1)-0.5 numel(temp_list)-1+1 numel(temp_list)-1+1], 'EdgeColor', 'r','LineWidth',2);
end

end

function if_plot_trial_trial_image(image_Z, metric)
imagesc(image_Z);
%axis image;
title(sprintf('tr-tr similarity ward sorted, %s', metric));
%caxis([-.2 .8]);
axis tight;
axis equal;
xlabel('Trials');
ylabel('Trials');
colorbar;

end