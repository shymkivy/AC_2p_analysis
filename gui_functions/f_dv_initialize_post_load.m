function f_dv_initialize_post_load(app)

app.DatasetDropDown.Items = app.data.dset_name_full;

app.trialtypeDropDown.Items = [{'all'; 'Freqs'; 'Context'; 'Context_flip';...
                    'Context_both'; 'Cont_dev'; 'Cont_dev_flip';...
                    'Ongoing'}; app.ops.context_types_labels];

max_planes = max(app.data.num_planes);
num_dsets = size(app.data,1);

for n_var = 1:numel(app.gui_ops.save_var_list_pl)
    var1 = app.gui_ops.save_var_list_pl{n_var};
    if ~sum(strcmpi(app.data.Properties.VariableNames, var1))
        app.data.(var1) = cell(num_dsets,max_planes);
    end
end

for n_var = 1:numel(app.gui_ops.save_var_list)
    var1 = app.gui_ops.save_var_list{n_var};
    if ~sum(strcmpi(app.data.Properties.VariableNames, var1))
        app.data.(var1) = cell(num_dsets,1);
    end
end

app.data.registered_data = cell(num_dsets,max_planes);

f_dv_load_dset_from_data(app);

end