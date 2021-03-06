trial_types = data.A1.trial_types_pr{1}(1:800);
traces = data.A1.tuning_all{1}.peak_tuning_full_resp.fr_peak_mag(:,1:800);

response = trial_types == 170;
predictors = traces';


SVMModel = fitcsvm(...
    predictors, ...
    response, ...
    'KernelFunction', 'linear', ...     % 'KernelFunction', 'polynomial'
    'PolynomialOrder', [], ...               % 2
    'KernelScale', 'auto', ...
    'BoxConstraint', 1, ...
    'Standardize', true, ...
    'ClassNames', [0; 1]);

% sv = SVMModel.SupportVectors;
% figure
% gscatter(X(:,1),X(:,2),y)
% hold on
% plot(sv(:,1),sv(:,2),'ko','MarkerSize',10)
% legend('versicolor','virginica','Support Vector')
% hold off


% 
% % linear or quadratic
% classificationDiscriminant = fitcdiscr(...
%     predictors, ...
%     response, ...
%     'DiscrimType', 'linear', ...        % 'DiscrimType', 'quadratic'
%     'Gamma', 0, ...
%     'FillCoeffs', 'off', ...
%     'ClassNames', [0; 1]);
% 
% % tree
% classificationTree = fitctree(...
%     predictors, ...
%     response, ...
%     'SplitCriterion', 'gdi', ...
%     'MaxNumSplits', 20, ...     % fine = 100; medium = 20; coarse = 4
%     'Surrogate', 'off', ...
%     'ClassNames', [0; 1]);
% 
% % naive beyes
% classificationNaiveBayes = fitcnb(...
%         predictors, ...
%         response, ...
%         'Kernel', 'Gaussian', ...
%         'Support', 'Unbounded', ...
%         'DistributionNames', distributionNames, ...
%         'ClassNames', [0; 1]);