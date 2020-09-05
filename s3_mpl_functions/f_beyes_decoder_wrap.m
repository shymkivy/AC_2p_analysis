function accuracy = f_beyes_decoder_wrap(predictors, response, tt, params)
KernelFunction = f_get_param(params, 'KernelFunction', 'gaussian');
KernelScale = f_get_param(params, 'KernelScale', 7); % auto fine = 1.8 medium= 7.1 coarse = 28
PolynomialOrder = f_get_param(params, 'PolynomialOrder');
kFold = f_get_param(params, 'KFold', 5);


%% split data
num_trials = numel(response);

cv_groups = f_make_crossval_groups(num_trials, kFold);

kFold = sum(logical(cv_groups.num_test_trials));

accuracy1 = zeros(kFold,1);
for n_cv = 1:kFold
    test_gr = cv_groups.test_trial_bool(:,n_cv);
    
    beyes_model = f_beyes_decoder_train(predictors(~test_gr,:), response(~test_gr), tt);
    
    beyes_prediction = f_beyes_decoder_predict(beyes_model, predictors(test_gr,:), tt);
    
    accuracy1(n_cv) = sum(beyes_prediction == response(test_gr))/sum(test_gr);
end

accuracy = mean(accuracy1);

end