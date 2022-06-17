function f_dv_ensemble_extract(app)

disp('Extracting ensembles...');

params = f_dv_gather_params(app);
params.cdata = f_dv_get_cdata(app);


%%
ens_out = f_dv_ensemble_extract_core(params);

%%
ddata_idx = strcmpi(app.ddata.dset_name_full, app.data.dset_name_full);
app.data(ddata_idx,:).ensembles{1} = ens_out;
app.ddata.ensembles{1} = ens_out;

disp('Done');

end