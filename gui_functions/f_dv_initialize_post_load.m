function f_dv_initialize_post_load(app)

app.DatasetDropDown.Items = app.data.experiment;

app.trialtypeDropDown.Items = [{'all'}; {'Freqs'}; {'Context'}; {'Ongoing'}; app.ops.context_types_labels];

max_planes = max(app.data.num_planes);

for n_var = 1:numel(app.gui_ops.save_var_list)
    var1 = app.gui_ops.save_var_list{n_var};
    if ~sum(strcmpi(app.data.Properties.VariableNames, var1))
        app.data.(var1) = cell(size(app.data,1),max_planes);
    end
end

app.data.registered_data = cell(size(app.data,1),max_planes);

f_dv_load_dset_from_data(app);
disp('Loaded')

end