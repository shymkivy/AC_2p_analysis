close all;
clear;

data_path = 'F:\AC_data\';
gui_dir = 'C:\Users\ys2605\Desktop\stuff\AC_2p_analysis';

experiment_type = 'tone_mmn'; % tone_mmn', 'FG_mmn', 'echo'

%%
addpath([gui_dir, '\analysis_functions\']);
addpath([gui_dir, '\gui_functions\']);
addpath([gui_dir, '\general_functions\']);
addpath([gui_dir, '\s3_mpl_functions\']);

ops = f_dset_ops(data_path);

ops.num_dsets_load = 20;

ops.experiment_type = experiment_type;

[data, ops, reg_struct] = f_load_data(ops);
params = ops.params;
params.paradigm = data.paradigm{1};
params.planes = 1;
params.n_pl = 1;
%%
params.region = 'All';      % all, all comb, a1, a2, uf, aaf
params.data_selection = 'All';   % all, mouse, dataset, plane
params.current_dset_idx = 1;
params.trial_window = [-0.05, 0.95];
params.trial_type = 'Context_both_comb';


params.convert_to_z = 1;
params.use_reg_data_labels = 1;
params.responsive_cells_type = 'peaks';
params.responsive_cells_select = 'resp marg';
params.responsive_thresh = 1;
params.stats_between = 'subdset';
params.stim_window = 'onset';
params.max_y_lim = 0;
params.min_y_lim = 0;
params.plot_stim = 1;
params.stim_freq_color = 4;
params.stim_transparancy = 0.2;
params.plot_super_deets = 0;

f_dv_plot_mmn(data, params, ops)

%%
params.region = 'All comb'; 
params.data_selection = 'All';
params.trial_type = 'Freqs -2';
params.trial_num_selection = 'min'; % all, median, mean, min
params.decoder_type = 'svm';        % svm, bayes, tree

[data_all, tt_all, reg_id, group_id] = f_dv_decoder_gather_data(data, params, ops);

dec_data = f_decoder_binwise_onevall(data_all, tt_all, params);

%%
params.region = 'All comb'; 
f_dv_decoder_onevall(data, params, ops);

%%
params.distance_method = 'cosine';  % cosine, euclidean, correlation, hamming, jaccard
params.do_similarity = 1;
params.plot_feature = 'peak resp mag z';
params.mat_tri = 'Ltri';    % Ltri - lower triangular, Utri - upper triangular, Full
params.colormap = 'gray';
f_dv_similarity_onevone(data, params, ops)

%% trial-to-trial analysis of CDR
params.distance_reference = 'pairwise';  % pairwise, zero, trial ave
f_dv_trial_to_trial_corr(data, params, ops)

%% needs ens analysis for this
params.sort_trials = 1;
params.sort_with_full_firing_rate = 1;
params.resort_by_ens = 1;
ddata = data(1,:);
cdata = f_dv_compute_cdata_mpl(ddata, params);
f_dv_trial_to_trial_w_full_rates(ddata, cdata, params, ops)


