%%%%%%%%%%%%%%%%%%%%%
%
%   V2 used a new struct format of saved data in previous step
%   
%   % editing deviance plots
%   % need to make spike inference better
%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear;
close all;

ops.processed_data_tag = 'processed_data';
ops.OA_output_tag = 'results_sort';
addpath([pwd '\s3_mpl_functions']);
addpath([pwd '\general_functions']);
%% List of files to load

%which files to analyze
ops.blah = 1;

%AC_files_preOA_7_31_19
%AC_files_OA_10_27_19
%AC_files_MPL_files_12_11_19
%AC_files_OA_3_16_20;
AC_files_OA_5_23_20

%% preprocessing parameters
% ------- Load params -----------
ops.normalize_firing_rate = 1;
ops.extra_SNR_thresh = 0; % 0 = no thresh
ops.redundent_to_analyze = 3;
ops.redundent_pool_trials = 2:8;
ops.dev_cells_ctx = 'ctx_tuned';      % options: 'all', 'ctx_tuned', 'tuned_all'
ops.remove_early_dev = 1;

% A1 A2 AAF DM/DA
ops.regions_to_analyze = {'A1','A2','DF','AAF'}; %,, ,,     % choose from fieldnames above 
% reg, flip, combined
ops.flip_to_analyze = [1 2 3];

% which type of infered signal you want to use
ops.signal_inference = 'MCMC'; % options: 'smooth_dfdt', 'MCMC', 'c_foopsi', 'df_f', 'raw'
% create another trace with extra smoothing for population anazysis
ops.signal_extra_smooth_sig = .100;   % in sec
ops.signal_extra_smooth_plot_examples = 0;

% computing threshholds
% the distribution of trial averaged data looks like cut off normal dist,
% so using median is better approximation for z score computation but it
% requires sampling data so it is slow. 
ops.stat.thresh_method = 'zscore_around_mean'; % options: 'ecdf_percentile', 'zscore_around_mean', 'zscore_around_median'
ops.stat.ecdf_percentile_thresh = 99;      % 95% 99.7
ops.stat.num_samples_drawn = 1000;
ops.stat.z_scores_thresh = 3;
ops.stat.z_scores_average_thresh = 1;
ops.stat.trials_to_sample = [1:10, 19, 20, 29, 30];%[1:18, 20, 21:28, 30]; % 1:10
ops.stat.reliability_thresh = .15;          % threshold for how many individual responses are required for being responsive
ops.stat.plot_examples = 0; % how many random examples to plot
ops.stat.save_est_samp = 'stat_save_4_9_20';
ops.stat.min_samp_size_z_normalization = 0;

% remove locomotion trials
ops.remove_loco_trials = 0;

% ---------------------Analysis window parameters--------------------
% define windows
ops.trial_window = [-.05, .95];       % in sec
ops.onset_window = [0.05, 0.5];
ops.offset_window = [0.55, 0.95];
ops.resp_window_time = [.05 .95];
ops.trial_window_long = [-5, 6];

% trial averaging baseline removal
ops.baseline_removal_trial_ave = 'kde'; % options: 'kde', 'mean', 'min', 'none' (default kde)

% % ----------------------------stim parameters-----------------------------
ops.stim.stim_isi = 0.5;
ops.stim.stim_duration = 0.5;

ops.stim.num_freqs = 10;
ops.stim.start_freq = 2000;
ops.stim.increase_factor = 1.5;


%%
ops = f_mpl_process_ops(ops);

%%
[data, ops] = f_mpl_load_data(ops);
%[data, ops] = f_mpl_load_data2(ops);

%%
data = f_mpl_preprocess_data(data, ops);

%% compute tunning
data = f_mpl_compute_tuning(data, ops);

%%
data = f_process_tuned_data(data, ops);

%%
disp('Done preprocesing');