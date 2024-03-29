function data_dim_est = f_dv_estimate_dim_pca_core(params)
%%
est_params_pca = params.est_params_pca;
%%
% tn_all = f_dv_get_trial_number(app);
% tt_all = app.ops.context_types_all(tn_all)';
% 
% stim_times = app.ddata.stim_frame_index{n_pl};
% mmn_freq = app.ddata.MMN_freq{1};
% trial_types = app.ddata.trial_types{1};

%%
firing_rate = cat(1,params.cdata.S_sm);
active_cells = sum(firing_rate,2) ~= 0;
firing_rate(~active_cells,:) = [];
%firing_rate = f_dv_get_firing_rate(params.cdata);

num_cells = size(firing_rate,1);
firing_rate = firing_rate(randperm(num_cells),:);

data_dim_est = f_ensemble_comp_data_dim2(firing_rate, est_params_pca);

end