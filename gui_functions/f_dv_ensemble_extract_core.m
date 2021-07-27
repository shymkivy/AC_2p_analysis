function ens_params = f_dv_ensemble_extract_core(app, params)

corr_dim = round(params.data_dim_pca.dimensionality_corr);

%% input paramseters for ensemble analysis
params2 = f_dv_ensemble_params(app, corr_dim);
ens_params = params2.ens_params;

%%
ens_params.vol_period = params.cdata.volume_period;

%%
firing_rate = params.cdata.S;

active_cells = sum(firing_rate,2) ~= 0;
firing_rate(~active_cells,:) = [];

%num_cells = size(firing_rate,1);
%firing_rate = firing_rate(randperm(num_cells),:);

%% extract ensambles
firing_rate_sm = f_smooth_gauss(firing_rate, ens_params.smooth_SD/ens_params.vol_period);
firing_rate_sm_norm = f_smooth_gauss(firing_rate_sm, ens_params.smooth_SD/ens_params.vol_period);
ens_out = f_ensemble_analysis_YS_raster(firing_rate_sm_norm, ens_params);

%% evaluate components
firing_rate_norm = f_normalize(firing_rate, ens_params.normalize);
acc_out_d = f_evaluate_ens_cv(ens_out, firing_rate_norm, ens_params);
ens_params.acc_out_d = acc_out_d;
ens_params.ens_out = ens_out;

end