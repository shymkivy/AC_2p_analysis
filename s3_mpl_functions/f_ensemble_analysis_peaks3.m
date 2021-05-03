function ens_out = f_ensemble_analysis_peaks3(firing_rate, params)
% input either 3D tials data (Cell x Time X Trial)
%           or 2D trial data (Cell x Trial)
%% parameters
%disp('Ensemble detection');
if ~exist('params', 'var') || isempty(params)
    params = struct;
end
%cond_name = f_get_param(params, 'cond_name', 'none');
%n_dset = f_get_param(params, 'n_dset', 0);
num_comps = f_get_param(params, 'num_comps', []);   % compute with dim est
normalize1 = f_get_param(params, 'normalize', 'norm_mean_std');  % 'norm_mean_std', 'norm_mean' 'none'
total_dim_thresh = f_get_param(params, 'total_dim_thresh', .7);
shuffle_method = f_get_param(params, 'shuffle_method', 'scramble');     % 'circ_shift' or 'scramble'
%corr_comp_thresh = f_get_param(params, 'corr_comp_thresh', .90);
use_LR_proj = f_get_param(params, 'use_LR_proj', 0);
ensamble_method = f_get_param(params, 'ensamble_method', 'nmf'); % 'svd', 'ICA', 'NMF', 'SPCA', 'tca', 'fa', 'gpfa'
ensamble_extraction = f_get_param(params, 'ensamble_extraction', 'thresh'); % clust 'thresh'
% NMF is best with thresh extration or clustering may be ok too
% ICA or SVD extraction best with clustering
plot_stuff = f_get_param(params, 'plot_stuff', 0);
% so far NMF and ICA seems to work
%plot_stuff_extra = 0;
%plot_ens_details = 0;


%%
ndims1 = ndims(firing_rate);
if ndims1 == 3
    [num_cells, ~, num_trials] = size(firing_rate);
    firing_rate = reshape(firing_rate, num_cells,[]);
elseif ndims1 == 2
    [~, num_trials] = size(firing_rate);
    %num_bins = 1;
end  

active_cells = sum(firing_rate,2) > 0;
firing_rate(~active_cells,:) = [];

firing_rate_norm = f_normalize(firing_rate, normalize1);

num_cells = size(firing_rate_norm,1);

%% dim reduction with SVD to calulate components number

% [U,S,V] = svd(firing_rate_norm);
% sing_val_sq = diag(S'*S);
% d_explained = sing_val_sq/sum(sing_val_sq)*100;
[d_coeff,d_score,~,~,d_explained,d_mu] = pca(firing_rate_norm');

%figure; plot(d_explained)
dimensionality_total_norm = sum(cumsum(d_explained)<(total_dim_thresh*100));
%[coeff,score,~,~,d_explained,~] = pca(firing_rate_norm');


%% repeat with not norm
[~,S2,~] = svd(firing_rate);
sing_val_sq2 = diag(S2'*S2);
d_explained2 = sing_val_sq2/sum(sing_val_sq2)*100;
%figure; plot(d_explained)
dimensionality_total = sum(cumsum(d_explained2)<(total_dim_thresh*100));

%% shuff and PCA

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
dimensionality_total_norm_shuff = mean(dim_total_shuff);
% eigenvalues below lower bound plus above upper should
% theoretically equal total number of neurons in all ensembles
%dimensionality_corr = sum(d_explained>prctile(max_lamb_shuff, corr_comp_thresh*100));
%dimensionality_corr = mean(sum(d_explained>max_lamb_shuff'));

comp_num_data = sum(d_explained>max_lamb_shuff');
dimensionality_corr = mean(sum(d_explained>max_lamb_shuff'))+std(comp_num_data);

if isempty(num_comps)
    num_comps = ceil(dimensionality_corr);
end
data_dim_est.dimensionality_total = dimensionality_total;
data_dim_est.dimensionality_first_comp_size = d_explained2(1);
data_dim_est.dimensionality_total_norm = dimensionality_total_norm;
data_dim_est.dimensionality_total_norm_shuff = dimensionality_total_norm_shuff;
data_dim_est.dimensionality_corr = dimensionality_corr;
data_dim_est.num_comps = num_comps;
data_dim_est.d_explained = d_explained(1:num_comps);
%data_dim_est.corr_comp_thresh = corr_comp_thresh;
data_dim_est.num_cells = num_cells;
data_dim_est.num_trials = num_trials;

%%
%SI_firing_rate = similarity_index(firing_rate_norm, firing_rate_norm);
%SI_firing_rate_shuff = similarity_index(firing_rate_shuff, firing_rate_shuff);

%firing_rate_LR = U(:,1:num_comps)*S(1:num_comps,1:num_comps)*V(:,1:num_comps)';
%SI_firing_rate_LR = similarity_index(firing_rate_LR, firing_rate_LR);
n_comp = 1:num_comps;
firing_rate_LR = (d_coeff(:,n_comp)*d_score(:,n_comp)'+d_mu')';

%% sort cells and trials
hc_params.method = params.hcluster_method;
hc_params.distance_metric = params.hcluster_distance_metric;
hc_params.plot_dist_mat = plot_stuff;
hc_params.plot_clusters = plot_stuff;
hc_params.num_clust = num_comps+1;
hc_params.title_tag = 'Coeffs (cells)';
hclust_out_cell = f_hcluster_wrap(d_coeff(:,n_comp), hc_params);
ord_cell = hclust_out_cell.dend_order;
ens_out.ord_cell = ord_cell;

sort_tr = 1;
if sort_tr
    d_score_norm = d_score(:,n_comp)./vecnorm(d_score(:,n_comp));
    hc_params.title_tag = 'Scores (trials)';
    hclust_out_tr = f_hcluster_wrap(d_score_norm, hc_params);
    ord_tr = hclust_out_tr.dend_order;
    ens_out.ord_tr = ord_tr;
end
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
        ens_out = f_ensemble_extract_clust(coeffs, scores, num_comps, [], params);
    elseif strcmpi(ensamble_extraction, 'thresh')
        [thresh_coeffs, thresh_scores] = f_ens_get_thresh(firing_rate_ensemb, coeffs, scores, num_comps, params);
        ens_out = f_ensemble_apply_thresh(coeffs, scores, thresh_coeffs, thresh_scores, num_comps);
    end

else
    ens_out.cells.clust_label = 0;
    ens_out.cells.ens_list = {(1:num_cells)'};
    ens_out.cells.clust_ident = zeros(num_cells,1);
    ens_out.trials.clust_label = 0;
    ens_out.trials.ens_list = {(1:num_trials)'};
    ens_out.trials.clust_ident = zeros(num_trials,1);
end

ens_out.cells.dend_order = ord_cell;
if sort_tr
    ens_out.trials.dend_order = ord_tr;
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
%     f_plot_comp_scatter(U(:,pl3_pc));
%     xlabel('pc 1');
%     ylabel('pc 2');
%     zlabel('pc 3');
%     title('U - cells');
% 
% 
%     %figure; plot3(V(:,1),V(:,2),V(:,3), 'o')
%     f_plot_comp_scatter(V(:,pl3_pc));
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
%         f_plot_comp_scatter(X(:,dims1));
%         xlabel(sprintf('comp %d', dims1(1)));
%         ylabel(sprintf('comp %d', dims1(2)));
%         zlabel(sprintf('comp %d', dims1(3)));
%         title(sprintf('%s, dset%d, trial type pcs scores(trials)',params.cond_name, params.n_dset));
%     end
% end

end

