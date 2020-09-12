function accuracy1 = f_svm_decoder(predictors, response, tt, params)
KernelFunction = f_get_param(params, 'KernelFunction', 'gaussian');
KernelScale = f_get_param(params, 'KernelScale', 7); % auto fine = 1.8 medium= 7.1 coarse = 28
PolynomialOrder = f_get_param(params, 'PolynomialOrder');
kFold = f_get_param(params, 'KFold', 5);

if ~strcmpi(KernelFunction, 'gaussian')
    KernelScale = [];
end

if ~strcmpi(KernelFunction, 'polynomial')
    PolynomialOrder = [];
end

%%
% SVMModel = fitcsvm(...
%     predictors(1:38,:), ...
%     response(1:38), ...
%     'KernelFunction', 'cosineKernel', ...     % 'KernelFunction', 'polynomial' 'gaussian'
%     'PolynomialOrder', [], ...               % 2
%     'KernelScale', [], ...
%     'BoxConstraint', 1, ...
%     'Standardize', true, ...
%     'ClassNames', tt(1:2)');

template = templateSVM(...
    'KernelFunction', KernelFunction, ...
    'PolynomialOrder', PolynomialOrder, ...
    'KernelScale', KernelScale, ...      % auto fine = 1.8 medium= 7.1 coarse = 28
    'BoxConstraint', 1, ...
    'Standardize', true);
SVMModel = fitcecoc(...
    predictors, ...
    response, ...
    'Learners', template, ...
    'Coding', 'onevsone', ...
    'ClassNames', tt');

%% cross validation 
model1 = SVMModel;

%x2 = predict(SVMModel, predictors);
%sum(x2 == response)/120

partitionedModel = crossval(model1, 'KFold', kFold);

%[validationPredictions, validationScores] = kfoldPredict(partitionedModel);
%sum(validationPredictions == response)/120


% Compute validation accuracy
accuracy1 = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');


end