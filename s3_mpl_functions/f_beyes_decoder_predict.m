function beyes_prediction = f_beyes_decoder_predict(beyes_model, test_data, tt)
vec_means = beyes_model.vec_means;
vec_dist_mean = beyes_model.vec_dist_mean;
vec_dist_std = beyes_model.vec_dist_std;

% decode identity based on the model
num_trials = size(test_data,1);
beyes_prediction = zeros(num_trials,1);

% test dist
test_dist = pdist2(vec_means, test_data, 'cosine');

for n_trial = 1:num_trials
    trial_p = exp(-((test_dist(:,n_trial) - vec_dist_mean).^2)./(2*vec_dist_std.^2));
    [~, dec_ident] = max(trial_p);
    beyes_prediction(n_trial) = tt(dec_ident);
end

end