function f_dv_ensemble_stats(app)

ddata = app.ddata;
ens_params = ddata.ensembles{1};

%%
if isempty(ens_params)
    disp('Run ensemble analysis first')
else
    
    cdata = f_dv_get_cdata(app);
    firing_rate = cat(1,cdata.S_sm);
    active_cells = sum(firing_rate,2) ~= 0;
    firing_rate(~active_cells,:) = [];
    
    %% shuffled data
    ens_params.acc_shuff_reps = 20;
    
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

    acc_out_full = cat(1,acc_out_s{:});
    thresh = prctile(acc_out_full(:), 1);
    
    
    ens_params.acc_out_shuff = acc_out_s;
    ens_params.acc_out_thresh = thresh;
    ens_params.accepted_ensembles = ens_params.acc_out_d<thresh;
    
    ens_stats = ens_params;
    
    ddata_idx = strcmpi(app.ddata.dset_name_full, app.data.dset_name_full);
    app.data(ddata_idx,:).ensemble_stats{1} = ens_stats;
    app.ddata.ensemble_stats{1} = ens_stats;

    app.numsigensEditField.Value = sum(ens_params.acc_out_d<thresh);
    
    disp('Done');

end
   
end