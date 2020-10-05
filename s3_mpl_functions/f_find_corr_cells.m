function corr_cells = f_find_corr_cells(raster_norm, dist_metric)
if ~exist('dist_metric', 'var') || isempty(dist_metric)
    dist_metric = 'correlation';    % correlation cosine
end
plot_stuff = 0;
corr_thresh = 95;
 
[num_cells, ~] = size(raster_norm);


dist_d = f_pdist_YS(raster_norm, dist_metric);
dist_d2 = dist_d + diag(ones(num_cells,1))*mean(squareform(dist_d));
pwcorr_d = mean(1-dist_d2,2);

num_reps = 100;
pwcorr_s_means = zeros(num_cells,num_reps);
for n_rep = 1:num_reps
    raster_norm_shuff = f_shuffle_data(raster_norm);
    dist_s = f_pdist_YS(raster_norm_shuff, dist_metric);
    dist_s2 = dist_s + diag(ones(num_cells,1))*mean(squareform(dist_s));
    pwcorr_s_means(:,n_rep) = mean(1-dist_s2);
end

dist_s_thresh_up = prctile(pwcorr_s_means, corr_thresh,2);
dist_s_thresh_down = prctile(pwcorr_s_means, 100-corr_thresh,2);

if plot_stuff
    figure; hold on; axis tight;
    plot(pwcorr_d);
    plot(dist_s_thresh_up, '--r');
    plot(dist_s_thresh_down, '--g');
    legend('mean pw corr', [num2str(corr_cell_thresh_percent) '% shuff thresh'], [num2str(100-corr_cell_thresh_percent) '% shuff thresh'])
    xlabel('cells'); ylabel('mean pw corr');
    title('correlated cells selection');
end

corr_cells = pwcorr_d>dist_s_thresh_up;

% 
% num_shuff_corr = 100;
% num_reps = 10;
% 
% for n_cell = 1
%     thresh_list = zeros(num_reps,2);
%     fprintf('Wait: ');
%     for n_rep = 1:num_reps
%         curr_cell = raster_norm(n_cell,:);
%         raster_rem = raster_norm;
%         raster_rem(n_cell,:) = [];
%         corr_out = 1 - f_pdist2_YS(curr_cell, raster_rem, dist_metric);
%         shuff_curr_cell = zeros(num_shuff_corr, num_cells-1);
%         for n_shuff = 1:num_shuff_corr
%             curr_cell_s = f_shuffle_data(curr_cell);
%             shuff_curr_cell(n_shuff,:) = 1 - f_pdist2_YS(curr_cell_s, raster_rem, dist_metric);
%         end
%         shuff_raster = zeros(num_shuff_corr, num_cells-1);
%         for n_shuff = 1:num_shuff_corr
%             raster_rem_s = f_shuffle_data(raster_rem);
%             shuff_raster(n_shuff,:) = 1 - f_pdist2_YS(curr_cell, raster_rem_s, dist_metric);
%         end
%         thresh_list(n_rep,1) = prctile(shuff_curr_cell(:), corr_thresh);
%         thresh_list(n_rep,2) = prctile(shuff_raster(:), corr_thresh);
%         fprintf('-%.0f%%-', n_rep/num_reps*100);
%     end
%     fprintf('\nDone\n');
% end
% 
% 
% sum(corr_out>mean(thresh_list(:,2)))
% 
% figure; hold on;
% histogram(corr_out,'Normalization','pdf')
% histogram(shuff_curr_cell(:),'Normalization','pdf')
% histogram(shuff_raster(:),'Normalization','pdf')
% 
% figure; hold on;
% histogram(thresh_list(:,1),'Normalization','pdf')
% histogram(thresh_list(:,2),'Normalization','pdf')
% 

end