function ens_out = f_ensemble_analysis_peaks3(trial_data, params, ops)
% input either 3D tials data (Cell x Time X Trial)
%           or 2D trial data (Cell x Trial)
%% parameters
disp('Ensemble detection');
cond_name = f_get_param(params, 'cond_name', 'none');
n_dset = f_get_param(params, 'n_dset', 0);
ens_list_gt = f_get_param(params, 'mark_cell_types');
trial_list_gt = f_get_param(params, 'mark_trial_types');
num_comps = f_get_param(params, 'num_comps', []);
normalize = f_get_param(params, 'normalize', 1);
total_dim_thresh = f_get_param(params, 'total_dim_thresh', .7);
shuffle_method = f_get_param(params, 'shuffle_method', 'scramble');     % 'circ_shift' or 'scramble'

var_thresh_prc = f_get_param(ops.ensemb, 'pca_var_thresh', .95);
ensamble_method = f_get_param(ops.ensemb, 'method', 'ica'); % 'svd', 'ICA', 'NMF', 'SPCA', 'tca', 'fa', 'gpfa'
cluster_method = f_get_param(params, 'cluster_method', 'hclust');    % 'hclust' or 'gmm'

plot_stuff = f_get_param(params, 'plot_stuff', 0);
% so far NMF and ICA seems to work
%plot_stuff_extra = 0;
%plot_ens_details = 0;

if ~isempty(ens_list_gt)
    compare_ground_truth = 1;
else
    compare_ground_truth = 0;
end

%%
ndims1 = ndims(trial_data);

if ndims1 == 3
    [num_cells, num_bins, num_trials] = size(trial_data);
    trial_data = reshape(trial_data, num_cells,[]);
elseif ndims1 == 2
    [~, num_trials] = size(trial_data);
    num_bins = 1;
end  

active_cells = sum(trial_data,2) > 0;
trial_data(~active_cells,:) = [];

if normalize
    firing_rate_norm = trial_data - mean(trial_data,2);
    firing_rate_norm = firing_rate_norm./std(firing_rate_norm,[],2); 
    %firing_rate_cont(isnan(firing_rate_cont)) = 0;
else
    firing_rate_norm = trial_data;
end

num_cells = size(firing_rate_norm,1);

%%
%firing_rate_norm_3d = reshape(firing_rate_norm, num_cells, num_bins, num_trials);

%% dim reduction with PCA to calulate components number

[U,S,V] = svd(firing_rate_norm);
sing_val_sq = diag(S'*S);
d_explained = sing_val_sq/sum(sing_val_sq)*100;
%figure; plot(d_explained)
dimensionality_total = sum(cumsum(d_explained)<(total_dim_thresh*100));

%% plotting
if plot_stuff
    for n_pc = 1:3
        im1 = (U(:,n_pc)*S(n_pc,n_pc)*V(:,n_pc)');
        figure;
        imagesc(im1);
        title(sprintf('comp %d', n_pc));
        if compare_ground_truth
            f_plot_cell_indicator(im1, ens_list_gt, ops);
            f_plot_trial_indicator2(im1, trial_list_gt, 1, ops);
        end
    end

    pl3_pc = 1:3;
    %figure; plot3(U(:,1),U(:,2),U(:,3), 'o')
    f_plot_comp_scatter(U(:,pl3_pc), ens_list_gt, ops);
    xlabel('pc 1');
    ylabel('pc 2');
    zlabel('pc 3');
    title('U - cells');


    %figure; plot3(V(:,1),V(:,2),V(:,3), 'o')
    f_plot_comp_scatter(V(:,pl3_pc), trial_list_gt, ops);
    xlabel('pc 1');
    ylabel('pc 2');
    zlabel('pc 3');
    title('V - trials')
end

% disp('fantastico no?');
% url= 'https://twitter.com/yusterafa';
% web(url)
%% shuff and PCA

num_reps = 20;
data_thresh = zeros(num_reps,1);
dim_total_shuff = zeros(num_reps,1);
for n_rep = 1:num_reps
    firing_rate_shuff = f_shuffle_data(firing_rate_norm, shuffle_method);
    [~,s_S,~] = svd(firing_rate_shuff);
    s_sing_val_sq = diag(s_S'*s_S);
    s_explained = s_sing_val_sq/sum(s_sing_val_sq)*100;
    
    dim_total_shuff(n_rep) = sum(cumsum(s_explained)<(total_dim_thresh*100));
    
    %[~,~,~,~,s_explained,~] = pca(firing_rate_shuff');
    data_thresh(n_rep) = prctile(s_explained, var_thresh_prc*100); % ss_explained or s_explained
end
dimensionality_total_shuff = mean(dim_total_shuff);
% eigenvalues below lower bound plus above upper should
% theoretically equal total number of neurons in all ensembles
dimensionality_corr = sum(sum(d_explained>data_thresh'))/num_reps;

if isempty(num_comps)
    num_comps = max(ceil(dimensionality_corr),1);
end
data_dim_est.dimensionality_total = dimensionality_total;
data_dim_est.dimensionality_total_shuff = dimensionality_total_shuff;
data_dim_est.dimensionality_corr = dimensionality_corr;
data_dim_est.num_comps = num_comps;
data_dim_est.d_explained = d_explained(1:num_comps);
data_dim_est.var_thresh_prc = var_thresh_prc;
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
%use_projection_mat = 1;
% if use_projection_mat
%     firing_rate_norm_dred = 
% end

%% real data 

num_LR_comps = num_comps;
if strcmpi(ensamble_method, 'nmf')
    num_LR_comps = round(num_comps*1.5);
end

[dred_factors1, ~] = f_dred_train2(firing_rate_LR, num_LR_comps, num_trials, ensamble_method);
[coeffs, scores] = f_dred_get_coeffs(dred_factors1);

%% Visualize traces


%% can either use the low rank trace or the scores (trial comps * singular vals)
%X = firing_rate_LR'
X = scores';

%% visualize
if plot_stuff
    num_plots = ceil(num_LR_comps/3);
    for n_plt = 1:num_plots
        dims1 = rem([1 2 3]+(n_plt-1)*3-1, num_LR_comps)+1;
        f_plot_comp_scatter(X(:,dims1), trial_list_gt, ops);
        xlabel(sprintf('pc %d', dims1(1)));
        ylabel(sprintf('pc %d', dims1(2)));
        zlabel(sprintf('pc %d', dims1(3)));
        title(sprintf('%s, dset%d, trial type pcs scores(trials)',params.cond_name, params.n_dset));
    end
end

%% measure similarity between sequential trials vs shuffled order
%f_measure_temp_similarity(scores, trial_types, tt_to_dred);

%% cluster trials
if strcmpi(cluster_method, 'hclust')
    %% cluster with hclust
    params2.method = 'cosine'; % cosine, ward
    params2.metric = 'cosine'; % cosine squaredeuclidean
    params2.plot_sm = plot_stuff;
    params2.num_clust = num_comps+1;
    clust_out = f_hcluster_trial3(firing_rate_LR', params2);
    %gscatter(X(:,1),X(:,2),hclust_out.clust_ident);
    
elseif strcmpi(cluster_method, 'gmm')
    %% clust with gmm, can find best regularizer
    optimize_reg_val = 0;
    if optimize_reg_val
        num_vals = 100;
        num_reps = 20;
        mean_acc = zeros(num_vals,num_reps);
        %rg_list = linspace(0.0001, 1, num_vals);
        rg_list = logspace(-5, 0, num_vals);
        for n_rg = 1:numel(rg_list)
            for n_rep = 1:num_reps
                params3.metric = 'cosine'; % cosine squaredeuclidean
                params3.RegularizationValue = rg_list(n_rg);
                params3.num_clust = num_comps+1;
                gmmclust_out = f_gmmcluster_trial(X, params3);
                %gscatter(X(:,1),X(:,2),gmmclust_out.clust_ident);
                eval_gmm = f_evaluate_ens_result(gmmclust_out.clust_ident, trial_list_gt, 0);
                mean_acc(n_rg, n_rep) = mean(eval_gmm.accuracy);
            end
        end
        f_rg = figure;
        shadedErrorBar(log10(rg_list), mean(mean_acc,2), std(mean_acc, [], 2)/sqrt(num_reps-1));
        xlabel('log10(rg val)')
        title('Click on optimal regression value');
        [x,~] = ginput(1);
        close(f_rg);
        params3.RegularizationValue = 10^x;
    else
        params3.RegularizationValue = 0.017;
    end
    params3.metric = 'cosine'; % cosine squaredeuclidean
    params3.num_clust = num_comps+1;
    clust_out = f_gmmcluster_trial(X, params3);
    [~, clust_out.dend_order] = sort(clust_out.clust_ident);
end
%%
clust_out2 = clust_out;
num_clust = clust_out2.num_clust;

if compare_ground_truth
    clust_eval_out = f_evaluate_ens_result(clust_out.clust_ident, trial_list_gt, plot_stuff);
    if plot_stuff
        suptitle(sprintf('%s sorting, mean acc = %.2f',cluster_method, mean(clust_eval_out.accuracy)));
    end
    clust_out2.clust_ident = clust_eval_out.ens_perm_ind(clust_out.clust_ident);
    fprintf('%s clustering with %.2f accuracy, %d clusters\n', cluster_method, mean(clust_eval_out.accuracy), num_comps+1);
end


%% reoder to make core ens first

clust_label = unique(clust_out2.clust_ident);
clust_centers = zeros(num_clust,num_comps);
clust_mag = zeros(num_clust,1);
for n_clust = 1:num_clust
    clust1 = clust_label(n_clust);
    clust_centers(n_clust,:) = mean(X(clust_out2.clust_ident == clust1,:));
    clust_mag(n_clust) = norm(clust_centers(n_clust,:));
end

[~, ens_order] = sort(clust_mag);

%%
ens_out = struct();
ens_out.clust_label = clust_label-1;
ens_out.num_clust = clust_out2.num_clust;
ens_out.dend_order_tr = clust_out2.dend_order';
ens_out.clust_ident_tr = ens_order(clust_out2.clust_ident)-1;
ens_out.clust_centers_tr = clust_centers(ens_order,:);
ens_out.clust_mag = clust_mag(ens_order);

ens_out.cell_num = zeros(num_clust,1);
for n_clust = 1:num_clust
    clust1 = ens_out.clust_label(n_clust);
    ens_out.cell_num(n_clust) = sum(ens_out.clust_ident_tr == clust1);
end


if plot_stuff
    f_plot_comp_scatter(X(:,1:3), ens_out.clust_ident_tr, ops);
    title(sprintf('Identified ensambles trials, %s', cluster_method));
end

%% get the ensemble cells

for n_clust = 2:num_clust
    cell_clust = zeros(num_cells,1);
    for n_cell = 1:num_cells
        [~, cell_clust(n_cell)] = min(sqrt(sum((coeffs(n_cell,:) - clust_centers).^2,2)));
    end
end







im1 = coeffs*clust_centers';



%pdist2(X,Y,Distance)


norm(coeffs(1,:) - clust_centers(2,:))

n1 = norm(coeffs - clust_centers(1,:))

figure; imagesc(coeffs*clust_centers')

figure; plot3(im1(:,1), im1(:,3), im1(:,4), 'o')

for n_clust = 1:(num_comps+1)
    clust1 = clust_label(n_clust);
    clust1_ind = clust_out.clust_ident == clust1;
    firing_rate_clust = firing_rate_norm*V(clust1_ind,1:num_comps)';%U(:,1:num_comps)*S(1:num_comps,1:num_comps)
end
find(clust1_ind)
find(ens_list_gt == clust_eval_out.ens_perm_ind(2))
im1 = firing_rate_norm(:,clust1_ind);
figure; imagesc(im1)
figure; imagesc(im1*im1')
figure; imagesc(firing_rate_LR_clust*firing_rate_LR_clust')

[U,S,V] = svd(im1)
figure; imagesc(V)


figure; plot(mean(im1,2))

