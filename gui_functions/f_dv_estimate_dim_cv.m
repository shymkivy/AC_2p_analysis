function f_dv_estimate_dim_cv(app)

if isempty(app.ddata.data_dim_pca{1})
    f_dv_estimate_dim_pca(app);
end

params = f_dv_gather_params(app);
params.ddata = app.ddata;
params.cdata = cat(1,params.ddata.cdata{:});

data_dim_cv = f_dv_estimate_dim_cv_core(params);

%cv_corr_dim = round(app.ddata.data_dim_pca{1}.dimensionality_corr);

%%
app.DimCVEditField.Value = data_dim_cv.dimensionality_corr;

max_planes = max(app.data.num_planes);
if ~sum(strcmpi(app.data.Properties.VariableNames, 'data_dim_cv'))
    app.data.data_dim_cv = cell(size(app.data,1),max_planes);
    app.ddata.data_dim_cv = cell(1,max_planes);
end

ddata_idx = strcmpi(app.ddata.dset_name_full, app.data.dset_name_full);
app.data(ddata_idx,:).data_dim_cv{1} = data_dim_cv;
app.ddata.data_dim_cv{1} = data_dim_cv;

f_dv_plot_dim_cv(app);

disp('Done')

end