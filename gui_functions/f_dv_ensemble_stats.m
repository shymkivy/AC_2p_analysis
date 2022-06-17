function f_dv_ensemble_stats(app)

ddata = app.ddata;

%%
if isempty(ddata.ensembles{1})
    disp('Run ensemble analysis first')
else
    
    params = f_dv_gather_params(app);
    params.ensembles = ddata.ensembles{1};
    params.cdata = f_dv_get_cdata(app);
    
    ens_stats = f_dv_ensemble_stats_core(params);
    
    ddata_idx = strcmpi(app.ddata.dset_name_full, app.data.dset_name_full);
    app.data(ddata_idx,:).ensemble_stats{1} = ens_stats;
    app.ddata.ensemble_stats{1} = ens_stats;

    app.numsigensEditField.Value = sum(ens_stats.accepted_ensembles);
    
    disp('Done');

end
   
end