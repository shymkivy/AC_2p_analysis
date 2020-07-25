function f_tsne(trial_peaks)


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

end