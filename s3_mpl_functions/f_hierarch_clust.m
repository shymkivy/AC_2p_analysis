function [dend_order, clust_ident] = f_hierarch_clust(data, dend_thresh)

if ~exist('dend_thresh', 'var')
    dend_thresh = 50;
end
method1 = 'ward';

figure;
[~, ~, dend_order] = dendrogram(linkage(data,method1), 1000,'ColorThreshold',dend_thresh,'Orientation','left');
figure;
Z = pdist(data(dend_order,:), 'cosine');
imagesc(1-squareform(Z));
axis image;
title('How many clusters?');
num_clust = input('How many clusters?');
title(['Cell to cell similarity sorted ' method1]);
%caxis([-.2 .8]);
axis equal tight;
xlabel('Cells');
ylabel('Cells');
colorbar;

figure;
[~, clust_ident, ~] = dendrogram(linkage(data,method1), num_clust,'ColorThreshold',1000,'Orientation','left');
close; 


ord1 = clust_ident(dend_order);
figure; hold on;
imagesc(1-squareform(Z));
axis image;
title(['Cell to cell similarity sorted ' method1]);
%caxis([-.2 .8]);
axis equal tight;
xlabel('Cells');
ylabel('Cells');
colorbar;
set(gca,'Ydir','reverse');
for n_clust = 1:num_clust
    temp_list = find(ord1 == n_clust);
    rectangle('Position',[temp_list(1) temp_list(1) numel(temp_list)-1 numel(temp_list)-1], 'EdgeColor', 'r','LineWidth',2)
end

end