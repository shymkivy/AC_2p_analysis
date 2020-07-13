%% run this after data is preprocessed with s31MPL

%% analysis params
% ----------------------- plot extra stuff --------------------------
ops.use_zscores = 1;


% for individual dsets
ops.plot.ctx_full_each_dset = 0;
ops.plot.dd_each_dset = 0;
ops.plot_combined = 1;

ops.tuning_plots = 0;
ops.ctx_plots = 0;
ops.spat_tuning_plot = 0;

ops.extra_clustering_plots = 0;
ops.centroid_plots = 0;

% under construction
ops.cluster_analysis = 0;
ops.population_analysis = 0;
ops.population_analysis_trials = 1;
ops.more_analysis = 0;

%% --------------------------- analysis params -------------------------
% ensemble analysis
ops.ensemb.method = 'tca'; %'tca'
ops.ensemb.pca_var_thresh = .95; % in percent
ops.ensemb.smooth_kernSD = 200;
ops.ensemb.select_upstates = false;
ops.ensemb.PCA_dim_reduction = false;
ops.ensemb.onset_time = .25;
ops.ensemb.offset_time = .75;

% dimenstionality reduction
ops.dred_params.trace_to_use = 'trials_specified'; % options: 'full' 'trials_all' 'trials_specified'
% cross validation method to compute dimensionality
ops.dred_params.do_cv = 0;
ops.dred_params.dim_estimate = 0;
ops.dred_params.randomize_trials = 1;
ops.dred_params.method_list = {'svd','nmf','tca'}; %  ,'gpfa','fa', 'spca'
ops.dred_params.num_comp = 1:5:50; %
ops.dred_params.kernSD = 200;    % ms
ops.dred_params.cv_num_folds = 4; % for train test
ops.dred_params.sort_trial_before = false;

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
f_mpl_plot_cond_details(data, ops);


%% ---------------------------population analysis--------------------------
if ops.population_analysis
    f_mpl_population_analysis(data, ops);
end


%% editing
if ops.population_analysis_trials
    f_mpl_population_analysis_trials(data, ops);
end
% 
%%
% ---------------------------clustering kmeans---------------------------
if ops.cluster_analysis
    f_mpl_cluster_analysis(data, ops);
end

f_mpl_cluster_analysis2(data, ops);

f_mpl_deviant_trials(data, ops);
% 
%     
% % ------------------------------more analysis----------------------------
% if ops.more_analysis
%     f_more_analysis(data, ops);
% end

%%
