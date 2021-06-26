function f_dv_compute_dim_red_pca(app)

disp('Estimating dimensionality...')

params = f_dv_gather_params(app);
params.cdata = app.cdata;

data_dim_est = f_dv_estimate_dim_pca_core(params);

app.DimpcaEditField.Value = data_dim_est.dimensionality_corr;

max_planes = max(app.data.num_planes);
if ~isfield(app.data, data_dim_est)
    app.data.data_dim_est = cell(size(app.data,1),max_planes);
    app.ddata.data_dim_est = cell(1,max_planes);
end

ddata_idx = strcmpi(app.ddata.experiment, app.data.experiment);
app.data(ddata_idx,:).data_dim_est{n_pl} = data_dim_est;
app.ddata.data_dim_est{n_pl} = data_dim_est;

disp('Done')

end