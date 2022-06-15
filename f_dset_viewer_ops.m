function f_dset_viewer_ops(app)

app.gui_ops.gui_save_dir = 'F:\AC_data\dset_viewer_save';

ops = struct();

%% preprocessing params

idx1 = 1;
ops.experiments(idx1).name = 'dream';
ops.experiments(idx1).data_path = 'F:\AC_data\caiman_data_dream3';
ops.experiments(idx1).save_mat_fname = 'dream_save_6_9_22.mat';
ops.experiments(idx1).save_reg_fname = '';

idx1 = 2;
ops.experiments(idx1).name = 'echo';
ops.experiments(idx1).data_path  = 'F:\AC_data\caiman_data_echo';
ops.experiments(idx1).save_mat_fname = '';
ops.experiments(idx1).save_reg_fname = '';

idx1 = 3;
ops.experiments(idx1).name = 'missmatch';
ops.experiments(idx1).data_path = 'F:\AC_data\caiman_data_missmatch';
ops.experiments(idx1).save_mat_fname = '';
ops.experiments(idx1).save_reg_fname = 'reg_save_6_10_21.mat';

idx1 = 4;
ops.experiments(idx1).name = 'missmatch_grating';
ops.experiments(idx1).data_path = '';
ops.experiments(idx1).save_mat_fname = '';
ops.experiments(idx1).save_reg_fname = '';

%%
% ----------- process ops params ----------
ops.regions_to_analyze = {'A1','AAF','A2','UF'}; %,, ,,     % choose from fieldnames above 

% ----------- load data params ----------
ops.waitbar = 1;

% ------- preprocess data params -----------
% ---------------------Analysis window parameters--------------------
% define windows
ops.plot_window = [-.5, 2];
ops.analysis_window = [-.05, .95];       % in sec
%ops.onset_window = [0.05, 0.5];
%ops.offset_window = [0.55, 0.95];
%ops.resp_window_time = [.05 .95];
% ops.trial_window_long = [-5, 6];

ops.normalize_firing_rate = 1;
ops.extra_SNR_thresh = 0; % 0 = no thresh


% ops.redundent_to_analyze = 3;
ops.redundent_pool_trials = 2:7;
% ops.dev_cells_ctx = 'ctx_tuned';      % options: 'all', 'ctx_tuned', 'tuned_all'
% ops.remove_early_dev = 1;

% 
% % A1 A2 AAF DF

% ops.flip_to_analyze = [1 2 3];  % 1 is regular, 2 is flip, 3 is combined
% 
% % which type of infered signal you want to use
% ops.signal_inference = 'MCMC'; % options: 'smooth_dfdt', 'MCMC', 'c_foopsi', 'df_f', 'raw'
% % create another trace with extra smoothing for population anazysis
% ops.signal_extra_smooth_sig = .100;   % in sec
% ops.signal_extra_smooth_plot_examples = 0;
% 
% % computing threshholds
% % the distribution of trial averaged data looks like cut off normal dist,
% % so using median is better approximation for z score computation but it
% % requires sampling data so it is slow.
% ops.stat.thresh_method = 'zscore_around_mean'; % options: 'ecdf_percentile', 'zscore_around_mean', 'zscore_around_median'
% ops.stat.ecdf_percentile_thresh = 99;      % 95% 99.7
% ops.stat.num_samples_drawn = 100;
% ops.stat.z_scores_thresh = 3;
% ops.stat.z_scores_average_thresh = 1; % average across different z thresh to control diff sample size
% ops.stat.trials_to_sample = [1:10, 19, 20, 29, 30];%[1:18, 20, 21:28, 30]; % 1:10
% ops.stat.reliability_thresh = .15;          % threshold for how many individual responses are required for being responsive
% ops.stat.min_samp_size_z_normalization = 0;
% ops.stat.plot_examples = 0; % how many random examples to plot
% % remove locomotion trials
% ops.remove_loco_trials = 0;
% 

% % trial averaging baseline removal
% ops.baseline_removal_trial_ave = 'kde'; % options: 'kde', 'mean', 'min', 'none' (default kde)
% 
% % % ----------------------------stim parameters-----------------------------
% ops.stim.stim_isi = 0.5;
% ops.stim.stim_duration = 0.5;
% 
% ops.stim.num_freqs = 10;
% ops.stim.start_freq = 2000;
% ops.stim.increase_factor = 1.5;

%% stat default params
stats.stat_method = 'shuff_pool'; % 'shuff_pool', 'shuff_locwise', 'z_thresh'
stats.stat_source = 'Freqs_dd'; % 'All', 'Freqs', 'Freqs_dd'
stats.z_thresh = 3;
stats.peak_bin_time = 0.250; % in sec
stats.num_shuff_samp = 2000;
stats.base_resp_win = [-1 3];
stats.loco_thresh = 99; % in percent;

%% est dim pca default params

est_params_pca.normalize = 'norm_mean_std'; % 'norm_mean_std', 'norm_mean' 'none'
est_params_pca.shuffle_method = 'circ_shift'; % 'circ_shift', 'scramble'
est_params_pca.dim_est_num_reps = 50;
est_params_pca.total_dim_thresh = 0.7;
est_params_pca.plot_stuff = 0;

%% dim est CV

est_params_cv.ensemble_method = 'pca';              % options: svd, pca (faster than svd), nmf, ica                % SVD is most optimal for encoding, NMF rotates components into something that is real and interpretable
est_params_cv.normalize = 'norm_mean_std'; % **'norm_mean_std'**, 'norm_mean' 'none'   % either way, need to normalize the power of signal in each cell, otherwise dimred will pull out individual cells
est_params_cv.shuffle_data_chunks = 1;   % 1 or 0, keeping cell correlations   % if the sequence of trial presentation contains information, you will need to shuffle. Also need to do in chunks because adjacent time bins are slightly correlated
% ---- input one or range of values to estimate across following

est_params_cv.smooth_SD_center = 0;       % larger window will capture 'sequences' of ensembles, if window is smaller than optimal, you will end up splitting those into more components
est_params_cv.smooth_SD_range = 0;
est_params_cv.smooth_SD_count = 1;

est_params_cv.num_comp_center_around_dim_pca = 1;
est_params_cv.num_comp_center = 8;     
est_params_cv.num_comp_range = 4;
est_params_cv.num_comp_count = 9;

est_params_cv.reps = 2;              % how many repeats per param 
est_params_cv.include_shuff_version = 0;

%% ensemble extract default params

ens_params.ensemble_method = 'nmf'; % options: svd, **nmf**, ica     % here NMF is
ens_params.num_comp_use_dim_pca = 1;
%ens_params.num_comp = corr_dim;
ens_params.num_comp = 10;
ens_params.smooth_SD = 0; % 110 is better?
ens_params.normalize = 'norm_mean_std'; % 'norm_mean_std', 'norm_mean' 'none'
ens_params.ensemble_extraction = 'thresh'; %  **'thresh'(only for nmf)** 'clust'(for all)
% --- for thresh detection (only nmf)
ens_params.ensemble_extraction_thresh = 'signal_z'; % 'shuff' 'signal_z' 'signal_clust_thresh'
ens_params.signal_z_thresh = 2;
ens_params.shuff_thresh_percent = 95;
% --- for clust detection and general sorting 
ens_params.hcluster_method = 'average';  % ward(inner square), **average**, single(shortest)     
ens_params.hcluster_distance_metric = 'cosine';  % none, euclidean, squaredeuclidean, **cosine**, hammilarity, rbf% for low component number better euclidean, otherwise use cosine
ens_params.corr_cell_thresh_percent = 95;   % to remove cells with no significant correlations
% --- other
ens_params.plot_stuff = 0;
ens_params.acc_shuff_reps = 20;

%%
gui_defaults.stat_method_options = {'shuff_pool', 'shuff_locwise', 'z_thresh'};
gui_defaults.stat_source_options = {'All', 'Freqs', 'Freqs_dd'};

gui_defaults.normalize_options = {'norm_mean_std', 'norm_mean' 'none'};
gui_defaults.shuffle_method_options = {'circ_shift', 'scramble'};

gui_defaults.ens_method_options = {'svd', 'nmf', 'ica'};
gui_defaults.ens_extraction_options = {'thresh', 'clust'};
gui_defaults.ens_extraction_thresh_options = {'shuff', 'signal_z', 'signal_clust_thresh'};

gui_defaults.est_dim_cv_method_options = {'pca', 'svd', 'nmf', 'ica'};

gui_defaults.hcluster_method_options = {'ward', 'average', 'single'};
gui_defaults.hcluster_distance_metric_options = {'none', 'euclidean', 'squaredeuclidean', 'cosine', 'hammilarity', 'rbf'};

%%
app.gui_ops.ops = ops;
app.gui_ops.ens_params = ens_params;
app.gui_ops.est_params_pca = est_params_pca;
app.gui_ops.est_params_cv = est_params_cv;
app.gui_ops.stats = stats;
app.gui_ops.gui_defaults = gui_defaults;
end