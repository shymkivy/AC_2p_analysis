function ens_out = f_ensemble_analysis_peaks3(trial_data, params, ops)
% input either 3D tials data (Cell x Time X Trial)
%           or 2D trial data (Cell x Trial)
%% parameters
disp('Ensemble detection');
if ~exist('params', 'var') || isempty(params)
    params = struct;
end
%cond_name = f_get_param(params, 'cond_name', 'none');
%n_dset = f_get_param(params, 'n_dset', 0);
num_comps = f_get_param(params, 'num_comps', []);   % compute with dim est
normalize = f_get_param(params, 'normalize', 'norm_full');  % 'norm_mean' 'norm_full'
total_dim_thresh = f_get_param(params, 'total_dim_thresh', .7);
shuffle_method = f_get_param(params, 'shuffle_method', 'scramble');     % 'circ_shift' or 'scramble'
corr_comp_thresh = f_get_param(params, 'corr_comp_thresh', .90);
use_LR_proj = f_get_param(params, 'use_LR_proj', 0);
ensamble_method = f_get_param(params, 'ensamble_method', 'nmf'); % 'svd', 'ICA', 'NMF', 'SPCA', 'tca', 'fa', 'gpfa'
ensamble_extraction = f_get_param(params, 'ensamble_extraction', 'clust'); % clust 'thresh'
% NMF is best with thresh extration or clustering may be ok too
% ICA or SVD extraction best with clustering
plot_stuff = f_get_param(params, 'plot_stuff', 0);
% so far NMF and ICA seems to work
%plot_stuff_extra = 0;
%plot_ens_details = 0;


%%
ndims1 = ndims(trial_data);

if ndims1 == 3
    [num_cells, ~, num_trials] = size(trial_data);
    trial_data = reshape(trial_data, num_cells,[]);
elseif ndims1 == 2
    [~, num_trials] = size(trial_data);
    %num_bins = 1;
end  

active_cells = sum(trial_data,2) > 0;
trial_data(~active_cells,:) = [];

if strcmpi(normalize, 'norm_full')
    firing_rate_norm = trial_data - mean(trial_data,2);
    firing_rate_norm = firing_rate_norm./std(firing_rate_norm,[],2); 
    %firing_rate_cont(isnan(firing_rate_cont)) = 0;
elseif strcmpi(normalize, 'norm_mean')
    firing_rate_norm = trial_data - mean(trial_data,2);
elseif strcmpi(normalize, 'none')
    firing_rate_norm = trial_data;
end

num_cells = size(firing_rate_norm,1);

%% dim reduction with SVD to calulate components number

[U,S,V] = svd(firing_rate_norm);
sing_val_sq = diag(S'*S);
d_explained = sing_val_sq/sum(sing_val_sq)*100;
%figure; plot(d_explained)
dimensionality_total = sum(cumsum(d_explained)<(total_dim_thresh*100));


%% shuff and PCA for num dim

num_reps = 200;
max_lamb_shuff = zeros(num_reps,1);
dim_total_shuff = zeros(num_reps,1);
for n_rep = 1:num_reps
    firing_rate_shuff = f_shuffle_data(firing_rate_norm, shuffle_method);
    [~,s_S,~] = svd(firing_rate_shuff);
    s_sing_val_sq = diag(s_S'*s_S);
    s_explained = s_sing_val_sq/sum(s_sing_val_sq)*100;
    
    dim_total_shuff(n_rep) = sum(cumsum(s_explained)<(total_dim_thresh*100));
    max_lamb_shuff(n_rep) = max(s_explained);
end
dimensionality_total_shuff = mean(dim_total_shuff);
% eigenvalues below lower bound plus above upper should
% theoretically equal total number of neurons in all ensembles
dimensionality_corr = sum(d_explained>prctile(max_lamb_shuff, corr_comp_thresh*100));
if isempty(num_comps)
    num_comps = ceil(dimensionality_corr);
end
data_dim_est.dimensionality_total = dimensionality_total;
data_dim_est.dimensionality_total_shuff = dimensionality_total_shuff;
data_dim_est.dimensionality_corr = dimensionality_corr;
data_dim_est.num_comps = num_comps;
data_dim_est.d_explained = d_explained(1:num_comps);
data_dim_est.corr_comp_thresh = corr_comp_thresh;
data_dim_est.num_cells = num_cells;
data_dim_est.num_trials = num_trials;

%%
%SI_firing_rate = similarity_index(firing_rate_norm, firing_rate_norm);
%SI_firing_rate_shuff = similarity_index(firing_rate_shuff, firing_rate_shuff);

firing_rate_LR = U(:,1:num_comps)*S(1:num_comps,1:num_comps)*V(:,1:num_comps)';
%SI_firing_rate_LR = similarity_index(firing_rate_LR, firing_rate_LR);

%%
%clust_est = f_hclust_estimate_num_clust(firing_rate_LR);

%%

if use_LR_proj
    firing_rate_ensemb = firing_rate_LR;
else
    firing_rate_ensemb = firing_rate_norm;
end

%% further reduce data 
if num_comps > 0
    num_LR_comps = num_comps;
    if strcmpi(ensamble_method, 'nmf')
        num_LR_comps = round(num_comps*1.5);
    end

    [dred_factors1, ~] = f_dred_train2(firing_rate_ensemb, num_LR_comps, ensamble_method, 0);
    [coeffs, scores] = f_dred_get_coeffs(dred_factors1);


    %%
    if strcmpi(ensamble_extraction, 'clust')
        ens_out = f_ensemble_extract_clust(coeffs, scores, num_comps, params, ops);
    elseif strcmpi(ensamble_extraction, 'thresh')
        [thresh_coeffs, thresh_scores] = f_ens_get_thresh(firing_rate_ensemb, coeffs, scores, num_comps, params);
        ens_out = f_ensemble_apply_thresh(coeffs, scores, thresh_coeffs, thresh_scores, num_comps);
    end

else
    ens_out.cells.clust_label = 0;
    ens_out.cells.ens_list = {(1:num_cells)'};
    ens_out.cells.clust_ident = zeros(num_cells,1);
    ens_out.cells.dend_order = f_hcluster(firing_rate_ensemb, 'cosine', 1);
    ens_out.trials.clust_label = 0;
    ens_out.trials.ens_list = {(1:num_trials)'};
    ens_out.trials.clust_ident = zeros(num_trials,1);
    ens_out.trials.dend_order = f_hcluster(firing_rate_ensemb', 'cosine', 1);
end


%ens_out.cells = f_get_clust_params(coeffs, ens_out.cells);
%ens_out.trials = f_get_clust_params(scores', ens_out.trials);
ens_out.data_dim_est = data_dim_est;

%% measure similarity between sequential trials vs shuffled order
%f_measure_temp_similarity(scores, trial_types, tt_to_dred);

%% Visualize traces

% if plot_stuff
%     for n_pc = 1:3
%         im1 = (U(:,n_pc)*S(n_pc,n_pc)*V(:,n_pc)');
%         figure;
%         imagesc(im1);
%         title(sprintf('comp %d', n_pc));
%     end
% 
%     pl3_pc = 1:3;
%     %figure; plot3(U(:,1),U(:,2),U(:,3), 'o')
%     f_plot_comp_scatter(U(:,pl3_pc), [], ops);
%     xlabel('pc 1');
%     ylabel('pc 2');
%     zlabel('pc 3');
%     title('U - cells');
% 
% 
%     %figure; plot3(V(:,1),V(:,2),V(:,3), 'o')
%     f_plot_comp_scatter(V(:,pl3_pc), [], ops);
%     xlabel('pc 1');
%     ylabel('pc 2');
%     zlabel('pc 3');
%     title('V - trials')
% end

% disp('fantastico no?');
% url= 'https://twitter.com/yusterafa';
% web(url)

% if plot_stuff
%     num_plots = ceil(num_LR_comps/3);
%     for n_plt = 1:num_plots
%         dims1 = rem([1 2 3]+(n_plt-1)*3-1, num_LR_comps)+1;
%         f_plot_comp_scatter(X(:,dims1), [], ops);
%         xlabel(sprintf('comp %d', dims1(1)));
%         ylabel(sprintf('comp %d', dims1(2)));
%         zlabel(sprintf('comp %d', dims1(3)));
%         title(sprintf('%s, dset%d, trial type pcs scores(trials)',params.cond_name, params.n_dset));
%     end
% end

end

