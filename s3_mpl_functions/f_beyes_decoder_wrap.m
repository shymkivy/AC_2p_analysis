function f_beyes_decoder_wrap(predictors, response, tt, temp_params(n_param_el))
KernelFunction = f_get_param(params, 'KernelFunction', 'gaussian');
KernelScale = f_get_param(params, 'KernelScale', 7); % auto fine = 1.8 medium= 7.1 coarse = 28
PolynomialOrder = f_get_param(params, 'PolynomialOrder');
kFold = f_get_param(params, 'KFold', 5);


%% split data






end