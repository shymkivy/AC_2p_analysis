function [out, dred_data] = f_dred_train(data, num_comp, method)
% data in data(cell, time, trial)
if ~exist('method', 'var') || isempty(method)
    method = 'SVD';
end

[num_cells, num_bins, num_trials] = size(data);


%% dim reduction

data_2d = reshape(data, num_cells, []);
data_means = zeros(num_cells,1);%mean(data_2d,2);
data_n2d = data_2d - data_means;


if strcmpi(method, 'svd')
    [U, S, V] = svd(data_n2d); 
    coeffs = U(:,1:num_comp);
    scores = S(1:num_comp,1:num_comp)*V(:,1:num_comp)'; 
    dred_data = coeffs * scores+data_means;
    dred_factors.scores = scores;
    dred_factors.coeffs = coeffs;
elseif strcmpi(method, 'nmf')
    data_means = zeros(num_cells,1);
    [d_W,d_H] = nnmf(data_2d,num_comp,'options',statset('Display','off'),'replicates',1,'algorithm','als');
    dred_data = d_W * d_H+data_means;
    dred_factors.d_W = d_W;
    dred_factors.d_H = d_H;
elseif strcmpi(method, 'fa')
    [faParams, ~] = fastfa(data_2d, num_comp, 'typ', 'fa');
    X = fastfa_estep(data_2d, faParams);
    dred_factors = faParams;
    % approximation of cov(x) x*x' ~ (Ph + L*L')
    % compute inverse (Ph + L*L')^-1 = iPh - iPh*L*(I + L'*iPh*L)^-1 * L'*iPh
    % dim reduced data = beta*x; beta = L'*(Ph + L*L')^-1
    %I = eye(num_comp);
    L = faParams.L;                        % (yDim-1) x 1
    %iPh = diag(1./faParams.Ph);
    %iPhL = iPh*L;   
    %MM = iPh - iPhL/(I+L'*iPhL)*iPhL';      % inverse of covariance appx;
    %beta = L'* MM;                          % beta is L'*covariance appx 
    %dred_data = L*beta*(data_s2d-faParams.d)+data_means;
    dred_data = L*X.mean+data_means;
elseif strcmpi(method, 'tca')
    data_sn3d = reshape(data_n2d, num_cells, num_bins, num_trials);
    est_factors = cp_als(tensor(data_sn3d),num_comp);
    dred_factors.t_factors = est_factors;
    dred_data = double(full(est_factors)) + data_means;
elseif strcmpi(method, 'gpfa') 
    seq = struct();
    for n_tr = 1:num_trials
        seq(n_tr).trialId = n_tr;
        seq(n_tr).T = num_bins;
        seq(n_tr).segId = 1;
        seq(n_tr).y = squeeze(data(:,:,n_tr));
    end
    
    [faParams, ~] = fastfa(data_2d, num_comp, 'typ', 'fa');
    
    binWidth      = 20; % in msec
    startTau      = 100; % in msec
    startEps      = 1e-3;
    startParams.covType = 'rbf';
    
    % GP timescale
    % Assume binWidth is the time step size.
    startParams.gamma = (binWidth / startTau)^2 * ones(1, num_comp);
    % GP noise variance
    startParams.eps   = startEps * ones(1, num_comp);
    startParams.d = data_means;
    startParams.C = faParams.L;
    startParams.R = diag(faParams.Ph);
    % Define parameter constraints
    startParams.notes.learnKernelParams = true;
    startParams.notes.learnGPNoise      = false;
    startParams.notes.RforceDiagonal    = true;
    currentParams = startParams;
    
    [estParams, seqTrainCut, LLcut, iterTime] = em(currentParams, seq);
    %[seqTrain, LLtrain] = exactInferenceWithLL(seq, estParams);
    

    % Leave-neuron-out prediction on test data 
    if estParams.notes.RforceDiagonal
      seq = cosmoother_gpfa_viaOrth_fast(seq, estParams, 1:num_comp);
    else
      seq = cosmoother_gpfa_viaOrth(seq, estParams, 1:num_comp);
    end
    % Compute log-likelihood of test data
    %[blah, LLtest] = exactInferenceWithLL(seqTest, estParams);

    fn = sprintf('ycsOrth%.2d',num_comp);
    dred_data = ([seq.(fn)]);
    dred_factors.estParams = estParams;
end

%train_err = norm(data_2d(:) - dred_data(:))/norm(data_2d(:));
train_err = norm(data_2d(:) - dred_data(:))/numel(data_2d);
% if strcmpi(method, 'nmf')
%     err_norm = norm(data_s2d(:) - dimRdata(:))/norm(data_s2d(:));
% else
%     err_norm = norm(data__sn2d(:) - dimRdata(:))/norm(data__sn2d(:));
% end

dred_factors.means = data_means;
out.num_comp = num_comp;
out.method = method;
out.dred_factors = dred_factors;
out.train_err_sm = train_err;
end