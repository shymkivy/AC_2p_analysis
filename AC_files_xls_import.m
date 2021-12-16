%% load data speciefied in xlsx file from given directory 

% drive4 = 'F';
% drive5 = 'G';
% grive6 = 'E';

% ops.file_dir = 'E:\data\AC\AC_data_OA_3_16_20';
% AC_data = readtable('AC_data_list.xlsx');
% ops.paradigm_type = 'ammn'; 

ops.file_dir = 'C:\Users\ys2605\Desktop\stuff\AC_data\caiman_data';
AC_data = readtable('AC_data_list_echo.xlsx');
ops.paradigm_type = 'cont'; % 'ammn' 'freq_grating' 'cont'

%%

use_dset = AC_data.im_use_dset;
use_dset(isnan(use_dset)) = 0;

AC_data.mpl(isnan(AC_data.mpl)) = 0;

%use_dset(AC_data.mpl<2) = 0;
%use_dset(AC_data.mpl>1) = 0;

AC_data = AC_data(logical(use_dset),:);
AC_data = AC_data(strcmpi(AC_data.paradigm,ops.paradigm_type),:);

ops.file_names.A1 = AC_data.experiment(strcmpi(AC_data.area, 'A1'));
ops.file_names.AAF = AC_data.experiment(strcmpi(AC_data.area, 'AAF'));
ops.file_names.A2 = AC_data.experiment(strcmpi(AC_data.area, 'A2'));
ops.file_names.UF = AC_data.experiment(strcmpi(AC_data.area, 'UF'));

ops.AC_data = AC_data;
clear AC_data use_dset;


%% preprocessing parameters
% ------- Load params -----------
ops.normalize_firing_rate = 1;
ops.extra_SNR_thresh = 0; % 0 = no thresh
ops.redundent_to_analyze = 3;
ops.redundent_pool_trials = 2:7;
ops.dev_cells_ctx = 'ctx_tuned';      % options: 'all', 'ctx_tuned', 'tuned_all'
ops.remove_early_dev = 1;
ops.waitbar = 1;

% A1 A2 AAF DF
ops.regions_to_analyze = {'A1','A2','AAF','UF'}; %,, ,,     % choose from fieldnames above 
ops.flip_to_analyze = [1 2 3];  % 1 is regular, 2 is flip, 3 is combined

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
ops.stat.num_samples_drawn = 100;
ops.stat.z_scores_thresh = 3;
ops.stat.z_scores_average_thresh = 1; % average across different z thresh to control diff sample size
ops.stat.trials_to_sample = [1:10, 19, 20, 29, 30];%[1:18, 20, 21:28, 30]; % 1:10
ops.stat.reliability_thresh = .15;          % threshold for how many individual responses are required for being responsive
ops.stat.min_samp_size_z_normalization = 0;
ops.stat.plot_examples = 0; % how many random examples to plot
% remove locomotion trials
ops.remove_loco_trials = 0;

% ---------------------Analysis window parameters--------------------
% define windows
ops.trial_window = [-.05, 1];       % in sec
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
