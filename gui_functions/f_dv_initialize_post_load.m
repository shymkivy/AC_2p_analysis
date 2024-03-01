function f_dv_initialize_post_load(app)

app.DatasetDropDown.Items = app.data.dset_name_full;

app.trialtypeDropDown.Items = [{'All'; 'all_no_cont'; 'All_orig_sequnce'; 'Freqs'; 'Freqs -1'; 'Freqs -2'; 'Context'; 'Context_flip';...
                    'Context_both'; 'Context_both_comb'; 'Cont_dev'; 'Cont_dev_flip'; 'Cont_dev_both'; 'Cont_dev_both_comb';...
                    'Context +1'; 'Context_flip +1'; 'Context_both +1'; 'Context_both_comb +1';...
                    'Context +2'; 'Context_flip +2'; 'Context_both +2'; 'Context_both_comb +2';...
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

end