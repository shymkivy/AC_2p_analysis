function ens_stats = f_dv_ensemble_stats_core(params)

%% input paramseters for ensemble analysis
%ens_params = params.ens_params;
ensembles = params.ensembles;
ens_params = params.ens_params;
ens_params.quiet = 1;
vol_period = params.cdata.volume_period;

%%
firing_rate = cat(1,params.cdata.S_sm);
active_cells = sum(firing_rate,2) ~= 0;
firing_rate(~active_cells,:) = [];

%num_cells = size(firing_rate,1);
%firing_rate = firing_rate(randperm(num_cells),:);

%%
if ens_params.acc_shuff_reps
    ens_params_s = ens_params;
    ens_params_s.num_comp = 10;
    ens_params_s.plot_stuff = 0;
    acc_out_s = cell(ens_params.acc_shuff_reps,1);
    fprintf('Computing shuff ens n/%d', ens_params.acc_shuff_reps);
    for n_shuff = 1:ens_params.acc_shuff_reps
        firing_rate_s = f_shuffle_data(firing_rate);
        
        % estimation
        firing_rate_s_sm = f_smooth_gauss(firing_rate_s, ens_params.smooth_SD/vol_period);
        firing_rate_s_sm_norm = f_normalize(firing_rate_s_sm, ens_params_s.normalize);
        ens_out_s = f_ensemble_analysis_YS_raster(firing_rate_s_sm_norm, ens_params_s);
        
        % accuracy
        firing_rate_s_norm = f_normalize(firing_rate_s, ens_params_s.normalize);
        acc_out_s2 = f_evaluate_ens_cv(ens_out_s, firing_rate_s_norm, ens_params_s);
        acc_out_s{n_shuff} = acc_out_s2(acc_out_s2>0);
        fprintf('--%d', n_shuff);
    end
    fprintf('\nDone\n');
    
    acc_out_full = cat(1,acc_out_s{:});
    thresh = prctile(acc_out_full, 1);

    ens_stats.ens_params = ens_params;
    ens_stats.acc_out_shuff = acc_out_s;
    ens_stats.acc_out_thresh = thresh;
    ens_stats.accepted_ensembles = ensembles.acc_out_d<thresh;

end

end