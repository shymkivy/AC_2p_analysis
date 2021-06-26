function f_dv_initialize_post_load(app)

app.DatasetDropDown.Items = app.data.experiment;

app.trialtypeDropDown.Items = [{'all'}; {'Freqs'}; {'Context'}; app.ops.context_types_labels];

max_planes = max(app.data.num_planes);

if ~sum(strcmpi(app.data.Properties.VariableNames, 'stats'))
    app.data.stats = cell(size(app.data,1),max_planes);
end

if ~sum(strcmpi(app.data.Properties.VariableNames, 'data_dim_pca'))
    app.data.data_dim_pca = cell(size(app.data,1),max_planes);
end

if ~sum(strcmpi(app.data.Properties.VariableNames, 'data_dim_cv'))
    app.data.data_dim_cv = cell(size(app.data,1),max_planes);
end

if ~sum(strcmpi(app.data.Properties.VariableNames, 'ensembles'))
    app.data.ensembles = cell(size(app.data,1),max_planes);
end

f_dv_load_dset_from_data(app);
disp('Loaded')

end