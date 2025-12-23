close all;
clear;

data_path = 'F:\AC_data\';
pipeline_dir = 'C:\Users\ys2605\Desktop\stuff\AC_2p_analysis';

%%
addpath([pipeline_dir, '\analysis_functions\']); % addpath(genpath())
addpath([pipeline_dir, '\gui_functions\']);
addpath([pipeline_dir, '\general_functions\']);
%addpath([pipeline_dir, '\s3_mpl_functions\']);

ops = f_dset_ops(data_path);

ops.num_dsets_load = 10;
ops.experiment_type = 'tone_mmn'; % tone_mmn', 'FG_mmn', 'echo'

[data, ops, reg_struct] = f_load_data(ops);
params = ops.params;

%%
params.region = 'All';      % all, all comb, a1, a2, uf, aaf
params.data_selection = 'All';   % all, mouse, dataset, plane
params.trial_type = 'Context_both_comb';

params.trial_window = [-0.05, 0.95];

f_dv_plot_mmn(data, params, ops)

%%
params.region = 'All'; 
params.data_selection = 'All';
params.trial_type = 'Freqs -1';
params.trial_num_selection = 'min'; % all, median, mean, min
params.decoder_type = 'svm';        % svm, bayes, tree
params.smooth = false;
params.smooth_sigma = 150;

[firing_rates_trials, trial_types_all, plot_t, region_id, reg_labels, trial_group_id] = f_dv_decoder_gather_data(data, params, ops);

n_dset = 1;
figure();
imagesc(reshape(firing_rates_trials{n_dset},size(firing_rates_trials{n_dset},1), []));
ylabel("Neurons"); xlabel("Frames"); title(sprintf("Dataset %d", n_dset))

%%
dec_data = f_decoder_binwise_onevall(firing_rates_trials, trial_types_all, params);

f_plot_decoder_data(dec_data, plot_t);

f_plot_decoder_data_by_regions(dec_data, plot_t, region_id, reg_labels, trial_group_id, params, ops);

%%
params.region = 'All comb'; 
f_dv_decoder_onevall(data, params, ops);

%%
params.distance_method = 'cosine';  % cosine, euclidean, correlation, hamming, jaccard
params.do_similarity = 1;
params.plot_feature = 'peak resp mag z';
params.mat_tri = 'Ltri';    % Ltri - lower triangular, Utri - upper triangular, Full
params.colormap = 'gray';
params.planes = 1;
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


%%

%f_dv_ensless_single_trial_corr(app)


