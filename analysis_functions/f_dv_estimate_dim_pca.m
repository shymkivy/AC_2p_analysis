function f_dv_estimate_dim_pca(app)

disp('Estimating dimensionality...')

params = f_dv_gather_params(app);
params.cdata = f_dv_get_cdata(app);
dim_est_pca = f_dv_estimate_dim_pca_core(params);

app.DimpcaEditField.Value = dim_est_pca.dimensionality_corr;

ddata_idx = strcmpi(app.ddata.dset_name_full, app.data.dset_name_full);
app.data(ddata_idx,:).data_dim_pca = {dim_est_pca};
app.ddata.data_dim_pca = {dim_est_pca};

f_dv_update_params(app);
disp('Done')
end