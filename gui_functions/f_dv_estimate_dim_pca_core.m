function data_dim_est = f_dv_estimate_dim_pca_core(params)
%%
est_params.normalize = 'norm_mean_std';
est_params.dim_est_num_reps = 50;
est_params.plot_stuff = 0;

%%
% tn_all = f_dv_get_trial_number(app);
% tt_all = app.ops.context_types_all(tn_all)';
% 
% stim_times = app.ddata.stim_frame_index{n_pl};
% mmn_freq = app.ddata.MMN_freq{1};
% trig_window = app.working_ops.trial_num_baseline_resp_frames;
% trial_types = app.ddata.trial_types{1};

%%
firing_rate = f_dv_get_firing_rate(params.cdata);

data_dim_est = f_ensemble_comp_data_dim2(firing_rate, est_params);

end