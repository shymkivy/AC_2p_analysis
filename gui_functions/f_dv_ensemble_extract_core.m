function ens_out = f_dv_ensemble_extract_core(params)

ens_params = params.ens_params;
vol_period = params.cdata.volume_period;
ens_params.quiet = 1;

%%
firing_rate = cat(1,params.cdata.S_sm);
active_cells = sum(firing_rate,2) ~= 0;
firing_rate(~active_cells,:) = [];

%num_cells = size(firing_rate,1);
%firing_rate = firing_rate(randperm(num_cells),:);

%% extract ensembles
firing_rate_sm = f_smooth_gauss(firing_rate, ens_params.smooth_SD/vol_period);
firing_rate_sm_norm = f_smooth_gauss(firing_rate_sm, ens_params.smooth_SD/vol_period);
ens_out = f_ensemble_analysis_YS_raster(firing_rate_sm_norm, ens_params);

%% evaluate components
firing_rate_norm = f_normalize(firing_rate, ens_params.normalize);
acc_out_d = f_evaluate_ens_cv(ens_out, firing_rate_norm, ens_params);
ens_out.acc_out_d = acc_out_d;
ens_out.ens_params = ens_params;

end