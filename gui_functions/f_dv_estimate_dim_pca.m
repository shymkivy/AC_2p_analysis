function f_dv_estimate_dim_pca(app)

disp('Estimating dimensionality...')

params = f_dv_gather_params(app);
params.cdata = app.cdata;

dim_est_pca = f_dv_estimate_dim_pca_core(params);

app.DimpcaEditField.Value = dim_est_pca.dimensionality_corr;

max_planes = max(app.data.num_planes);
if ~sum(strcmpi(app.data.Properties.VariableNames, 'dim_est_pca'))
    app.data.dim_est_pca = cell(size(app.data,1),max_planes);
    app.ddata.dim_est_pca = cell(1,max_planes);
end

ddata_idx = strcmpi(app.ddata.experiment, app.data.experiment);
app.data(ddata_idx,:).dim_est_pca{n_pl} = dim_est_pca;
app.ddata.data_dim_est{n_pl} = dim_est_pca;

disp('Done')
end