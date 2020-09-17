function bayes_prediction = f_bayes_decoder_predict(bayes_model, test_data, tt)
vec_means = bayes_model.vec_means;
vec_dist_mean = bayes_model.vec_dist_mean;
vec_dist_std = bayes_model.vec_dist_std;

% decode identity based on the model
num_trials = size(test_data,1);
bayes_prediction = zeros(num_trials,1);

% test dist
test_dist = pdist2(vec_means, test_data, 'cosine');

for n_trial = 1:num_trials
    trial_p = exp(-((test_dist(:,n_trial) - vec_dist_mean).^2)./(2*vec_dist_std.^2));
    [~, dec_ident] = max(trial_p);
    bayes_prediction(n_trial) = tt(dec_ident);
end

end