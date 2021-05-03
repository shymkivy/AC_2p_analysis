%%%%%%%%%%%%%%%%%%%%%
%
%   AC MMN analysis pipeline step 3.2
%   
%       step 1 processing voltage data
%       step 2 is obsolete
%       step 3.1 loads data and precomputes things for 3.2
%       step 3.2 analyis
%%%%%%%%%%%%%%%%%%%%%%%%%

%% run this after data is preprocessed with s31MPL
close all

%% analysis params
% ----------------------- plot extra stuff --------------------------
ops.use_zscores = 1;

% for individual dsets
ops.plot.ctx_full_each_dset = 1;
ops.plot.dd_each_dset = 1;
ops.plot_combined = 1;
ops.plot.spatial_tuning_dset = 1;
ops.plot.tuning_dset = 1;

ops.tuning_plots = 0;
ops.ctx_plots = 0;
%% --------------------------- analysis params -------------------------
% ensemble analysis
ops.ensemb.method = 'nmf'; %'tca'
ops.ensemb.corr_comp_thresh = .90; % in percent
ops.ensemb.total_dim_thresh = .70;
ops.ensemb.normalize = 'norm_mean_std'; % 'norm_mean_std', 'norm_mean' 'none'

ops.ensemb.smooth_kernSD = 200;
ops.ensemb.select_upstates = false;
ops.ensemb.PCA_dim_reduction = false;
ops.ensemb.onset_time = .25;
ops.ensemb.offset_time = .75;

% dimenstionality reduction
ops.dred_params.saved_data_dir = 's3_saved_dred_data';
ops.dred_params.run_idx = 2;    % for saving different runs

% select what trial types to use for trial-trial analysis as a cell list
% options: [numeric], 'all', 'cont_all'
% options: 'mmn1', 'mmn2', 'mmn12', 'dd1', 'dd2', 'dd12'
% options: 'red1', 'red2', 'red12', 'cont1', 'cont2', cont12'
ops.dred_params.trial_types_to_dred = {'dd1' 'dd2'}; % 
ops.dred_params.use_responsive_cells = 1;

%--------------
ops.dred_params.do_hclust = 1;
ops.dred_params.hclust.num_clust = {1, 1, 1};
ops.dred_params.hclust.method = 'average';  % ward(inner square), average, single(shortest)
ops.dred_params.hclust.plot_metric = 'cosine';  % none, euclidean, squaredeuclidean, cosine, hammilarity
ops.dred_params.hclust.sort_raster = 1;
ops.dred_params.hclust.plot_hclust_cells = 1;
ops.dred_params.hclust.plot_hclust_trials = 1;
%-----------------
ops.dred_params.do_cv = 0;
ops.dred_params.do_dim_estimate = 1;
ops.dred_params.randomize_trials = 1;
ops.dred_params.method_list = {'svd','nmf','tca'}; %  ,'gpfa','fa', 'spca'
ops.dred_params.num_comp = 1:2:20; %
ops.dred_params.kernSD = 200;    % ms
ops.dred_params.cv_num_folds = 4; % for train test
ops.dred_params.sort_trial_before = false;
%-----------------
ops.dred_params.do_ensamble_analysis = 1;
%
ops.dred_params.do_trial_dist = 1;
ops.dred_params.trial_types_for_dist = {[1:10], 'mmn12'};
% k-means param 

% how many clusters to use 
ops.num_clusters = [3 3 3];
% change order of clusters instead of descending by size, otherwise leave empty
ops.custom_clust_order = [2 1 3]; %[]; %

% normalization before dim reduction, yes or no
ops.normalize_before = 0;

% dimensionality reduction method
% 1 - take max, - this gets rid of time info
% 2 - PCA
% 3 - bin time
ops.dim_red_method = 1;

% for PCA dimensionality reduction
ops.num_scores_to_use = [3 2 2]; % C R D

% normalization after dim reduction
% 1 = none
% 2 = max(C,R,D)
% 3 = circle
% 4 = by std of each condition
ops.norm_after_method = 2;

%% Plot things
%f_mpl_plot_cell_details(data, ops);

%%
f_mpl_plot_dset_details(data, ops);

%%
%f_mpl_plot_cond_details(data, ops);

%%
%f_fov_registration(data, ops);
%% plot reliability of all cells
%f_plot_cond_reliability(data, ops);

%% make model of fake ensembles
%f_model_raster(data, ops);

%% ensemble analysis of full raster
%f_mpl_population_analysis_raster(data, ops);

%% ensemble analysis of trials
%f_mpl_population_analysis_trials2(data, ops);

f_mpl_population_analysis_trials3(data, ops);

%%
%f_mpl_trial_trial_dist(data, ops);

%%
%f_mpl_decoder(data, ops);

%% ---------------------------clustering trial aves---------------------------
f_mpl_trial_ave_clustering(data, ops);

%% -----------------------------more analysis----------------------------
%f_more_analysis(data, ops);


%%
