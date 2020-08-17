function [dend_order, clust_ident, Z] = f_hcluster(data, method, num_clust)

dend_thresh = 50;

% if ~exist('fig_hand', 'var')
%     fig_hand{1} = figure;
%     fig_hand{2} = figure;
% end

if strcmpi(method, 'ward')
    Z = linkage(data,'ward');
elseif strcmpi(method, 'cosine')
    Z = linkage(pdist(data,'cosine'),'average'); % weighted
end

f1 = figure;
[~, ~, dend_order] = dendrogram(Z, 3000,'ColorThreshold',dend_thresh,'Orientation','left');
close(f1);


if exist('num_clust', 'var')
    clust_ident = cluster(Z, 'MaxClust', num_clust);
else
    clust_ident = ones(size(data,1),1);
end



end

