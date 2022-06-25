function f_dv_ensemble_extract(app)

disp('Extracting ensembles...');

params = f_dv_gather_params(app);
params.cdata = f_dv_get_cdata(app);

if strcmpi(app.numenstofindmethodDropDown.Value, 'corr fraction')
    if ~isempty(app.ddata.data_dim_pca{1})
        params.ens_params.num_comp = ceil(app.ddata.data_dim_pca{1}.num_comps * app.numenstofindfracEditField.Value);
    else
        disp('Compute dim est pca to use corr fraction version; using value instead')
        params.ens_params.num_comp = app.numenstofindvalEditField.Value;
    end
elseif strcmpi(app.numenstofindmethodDropDown.Value, 'Value')
    params.ens_params.num_comp = app.numenstofindvalEditField.Value;
end

%%
ens_out = f_dv_ensemble_extract_core(params);

%%
ddata_idx = strcmpi(app.ddata.dset_name_full, app.data.dset_name_full);
app.data(ddata_idx,:).ensembles{1} = ens_out;
app.ddata.ensembles{1} = ens_out;

disp('Done');

end