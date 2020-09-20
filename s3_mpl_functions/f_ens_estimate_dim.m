function dim_corr = f_ens_estimate_dim(raster, num_shuff_reps)

% [d1,d2] = size(raster);
% d_min = min(d1,d2);

if ~exist('num_shuff_reps', 'var') || isempty(num_shuff_reps)
    num_shuff_reps = 50;
end

[~,~,~,~,d_explained,~] = pca(raster');

% [~,S,~] = svd(raster);
% S = S(1:d_min,1:d_min);
% sing_val_sq = diag(S).^2;
% d_explained = sing_val_sq/sum(sing_val_sq(:))*100;

max_lamb_shuff = zeros(num_shuff_reps,1);
for n_rep = 1:num_shuff_reps
    smooth_raster_shuff = f_shuffle_data(raster);
    
    [~,~,~,~,s_explained,~] = pca(smooth_raster_shuff');
    
%     [~,s_S,~] = svd(smooth_raster_shuff);
%     s_S = s_S(1:d_min,1:d_min);
%     s_sing_val_sq = diag(s_S).^2;
%     s_explained = s_sing_val_sq/sum(s_sing_val_sq(:))*100;

    max_lamb_shuff(n_rep) = max(s_explained);
end

comp_num_data = sum(d_explained>max_lamb_shuff');
dim_corr = mean(comp_num_data);

end