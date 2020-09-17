function bayes_model_out = f_bayes_decoder_train(predictors, response, tt)
% beyes decoder with cosine similarity kernel

num_cells = size(predictors,2);

% first get mean population vectors
vec_means = zeros(numel(tt),num_cells);
for n_tt = 1:numel(tt)
   tt_ind = response==tt(n_tt);
   vec_means(n_tt,:) = mean(predictors(tt_ind,:));
end

% compute distance from mean to each trial vector
vec_dist = pdist2(vec_means, predictors, 'cosine');

% compute mean distance and std for each trial type to create model
vec_dist_mean = zeros(numel(tt),1);
vec_dist_std = zeros(numel(tt),1);
for n_tt = 1:numel(tt)
    tt_ind = response==tt(n_tt);
    vec_dist_mean(n_tt) = mean(vec_dist(n_tt,tt_ind));
    vec_dist_std(n_tt) = std(vec_dist(n_tt,tt_ind));
end

bayes_model_out.vec_means = vec_means;
bayes_model_out.vec_dist_mean = vec_dist_mean;
bayes_model_out.vec_dist_std = vec_dist_std;

end

