function f_dv_estimate_dim_pca(app)

disp('Estimating dimensionality...')

params = f_dv_gather_params(app);

if strcmpi(app.SelectdatagroupButtonGroup.SelectedObject.Text, 'plane')
    params.cdata = app.cdata(n_pl,:);
else
    params.cdata = app.cdata;
end

dim_est_pca = f_dv_estimate_dim_pca_core(params);

app.DimpcaEditField.Value = dim_est_pca.dimensionality_corr;
% 
% max_planes = max(app.data.num_planes);
% if ~sum(strcmpi(app.data.Properties.VariableNames, 'dim_est_pca'))
%     app.data.dim_est_pca = cell(size(app.data,1),max_planes);
%     app.ddata.dim_est_pca = cell(1,max_planes);
% end

app.ddata.dim_est_pca = app.DimpcaEditField.Value;

ddata_idx = strcmpi(app.ddata.experiment, app.data.experiment);
app.data(ddata_idx,:).dim_est_pca = dim_est_pca;
app.ddata.data_dim_est = dim_est_pca;

disp('Done')
end