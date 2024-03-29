function accuracy = f_ens_estimate_corr_dim_cv(firing_rate, params)
% estimate dimensionality of correlations with a cross-validation method
if ~exist('params', 'var') || ~isstruct(params)
    params = struct;
end
kFold = f_get_param(params, 'KFold', 5);
shuffle_data_chunks = f_get_param(params, 'shuffle_data_chunks', 1);
num_comp = f_get_param(params, 'num_comp');
ensemble_method = f_get_param(params, 'ensemble_method', 'pca');
smooth_SD = f_get_param(params, 'smooth_SD', 100);
vol_period = f_get_param(params, 'vol_period', 1000/30);

renormalize = 1;

[num_cells, num_bins] = size(firing_rate);

cv_groups = f_make_crossval_groups(num_bins, kFold);

% get components to explain 25%
if isempty(num_comp)
    [~,~,~,~,explained,~] = pca(firing_rate);
    num_comp = sum(cumsum(explained) < 10);
end


if shuffle_data_chunks
    chunk_size = ceil(num_bins/kFold/100);
    chunk_idx = 0:(chunk_size):num_bins;
    num_chunks = numel(chunk_idx)-1;
    chunk_order = randperm(num_chunks);
    test_order1 = (1:num_bins)';
    test_order2 = reshape(test_order1(1:chunk_idx(end)), chunk_size, num_chunks);
    test_order3 = test_order2(:,chunk_order);
    test_orter4 = test_order3(:);
    if chunk_idx(end) < num_bins
        test_order = [test_orter4; test_order1((chunk_idx(end)+1):num_bins)];
    else
        test_order = test_orter4;
    end
else
    test_order = (1:num_bins)';
end

firing_rate_sm = f_smooth_gauss(firing_rate, smooth_SD/vol_period);

if renormalize
    firing_rate_sm = f_normalize(firing_rate_sm, 'norm_mean_std');
end

train_err = zeros(kFold,1);
train_err_sm = zeros(kFold,1);
test_err = zeros(kFold,1);
test_err_sm = zeros(kFold,1);
test_norm = zeros(kFold,1);
test_norm_sm = zeros(kFold,1);
for n_cv = 1:kFold
    test_gr = cv_groups.test_trial_bool(test_order,n_cv);
    
    yTrain = firing_rate(:,~test_gr);
    yTrain_sm = firing_rate_sm(:,~test_gr);
    [dred_factors, ydred_data] = f_dred_train2(yTrain_sm, num_comp, ensemble_method, 0);
    %train_err(n_cv) = norm(yTrain(:) - ydred_data(:))/norm(yTrain(:));
    %train_err_sm(n_cv) = norm(yTrain_sm(:) - ydred_data(:))/norm(yTrain_sm(:));
    train_err(n_cv) = norm(yTrain(:) - ydred_data(:))/numel(yTrain);
    train_err_sm(n_cv) = norm(yTrain_sm(:) - ydred_data(:))/numel(yTrain_sm);
    
    yTest = firing_rate(:,test_gr);
    yTest_sm = firing_rate_sm(:,test_gr);
    % reconstruct from not smooth data, but smooth is similar it seems
    Ycs = f_dred_test(yTest, dred_factors.dred_factors, ensemble_method); % _sm
    %test_err(n_cv) = norm(yTest(:) - Ycs(:))/norm(yTest(:));
    %test_err_sm(n_cv) = norm(yTest_sm(:) - Ycs(:))/norm(yTest_sm(:));
     
%     yTest_base = yTest - mean(yTest,2);
%     yTest_sm_base = yTest_sm - mean(yTest_sm,2);
%     Ycs_base = Ycs - mean(Ycs,2);
    
    %fac1 = sum(yTest_base(:) .* Ycs_base(:)/norm(Ycs_base(:)))/norm(Ycs_base(:));
    
    % normalized error, so cell number doesnt matter.. if data not
    % normalized after smoothing, error will not be flat
    test_err(n_cv) = norm(yTest(:) - Ycs(:))/norm(ones(num_cells,1));% norm(ones(size(yTest,1),1)); % norm(yTest(:));
    %fac(n_cv) = fac1;
    
    %fac1 = sum(yTest_sm_base(:) .* Ycs_base(:)/norm(Ycs_base(:)))/norm(Ycs_base(:));
    %test_err_sm(n_cv) = norm(yTest_sm(:) - Ycs(:))/numel(yTest_sm);
    test_err_sm(n_cv) = norm(yTest_sm(:) - Ycs(:))/norm(ones(num_cells,1));
    %fac_sm(n_cv) = fac1;
    
    test_norm(n_cv) = norm(yTest(:))/norm(ones(num_cells,1));
    
    test_norm_sm(n_cv) = norm(yTest_sm(:))/norm(ones(num_cells,1));
end


accuracy.train_err = nanmean(train_err);
accuracy.train_err_sm = nanmean(train_err_sm);
accuracy.test_err = nanmean(test_err);
accuracy.test_err_sm = nanmean(test_err_sm);
accuracy.test_norm = mean(test_norm);
accuracy.test_norm_sm = mean(test_norm_sm);
end