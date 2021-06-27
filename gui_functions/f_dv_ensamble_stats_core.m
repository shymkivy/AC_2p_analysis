function ens_stats = f_dv_ensamble_stats_core(~, params)

%% input paramseters for ensemble analysis
ens_params = params.ensembles;
ens_params.acc_shuff_reps = 20;

%%
firing_rate = params.cdata.S;

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
        firing_rate_s_sm = f_smooth_gauss(firing_rate_s, ens_params.smooth_SD/ens_params.vol_period);
        firing_rate_s_sm_norm = f_normalize(firing_rate_s_sm, ens_params_s.normalize);
        ens_out_s = f_ensemble_analysis_YS_raster(firing_rate_s_sm_norm, ens_params_s);
        
        % accuracy
        firing_rate_s_norm = f_normalize(firing_rate_s, ens_params_s.normalize);
        acc_out_s{n_shuff} = f_evaluate_ens_cv(ens_out_s, firing_rate_s_norm, ens_params_s);
        fprintf('--%d', n_shuff);
    end
    fprintf('\nDone\n');

    try
        acc_out_full = cat(1,acc_out_s{:});
    catch
        acc_out_full = [];
        for ii = 1:numel(acc_out_s)
            acc_out_full = [acc_out_full; acc_out_s{ii}(:)];
        end
        disp('error in acc extract');
    end
    thresh = prctile(acc_out_full(:), 1);

    ens_params.acc_out_shuff = acc_out_s;
    ens_params.acc_out_thresh = thresh;
    ens_params.accepted_ensembles = ens_params.acc_out_d<thresh;
end

ens_stats.ensembles = params.ensembles;
ens_stats.acc_shuff_reps = ens_params.acc_shuff_reps;
ens_stats.acc_out_shuff = ens_params.acc_out_shuff;
ens_stats.acc_out_thresh = ens_params.acc_out_thresh;
ens_stats.accepted_ensembles = ens_params.accepted_ensembles;
end