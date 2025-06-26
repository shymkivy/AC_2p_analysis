function f_dv_initialize_post_load(app)

app.regdatapathEditField.Value = app.ops.fpath_reg_data;
app.matdatapathEditField.Value = app.ops.fpath_mat_data;

app.DatasetDropDown.Items = app.data.dset_name_full;

app.trialtypeDropDown.Items = [{'All'; 'all_no_cont'; 'All_orig_sequnce'; 'Freqs'; 'Freqs -1'; 'Freqs -2'; 'Context'; 'Context_flip';...
                    'Context_both'; 'Context_both_comb'; 'Cont_dev'; 'Cont_dev_flip'; 'Cont_dev_both'; 'Cont_dev_both_comb';...
                    'Context +1'; 'Context_flip +1'; 'Context_both +1'; 'Context_both_comb +1';...
                    'Context +2'; 'Context_flip +2'; 'Context_both +2'; 'Context_both_comb +2';...
                    'Ongoing'}; app.ops.context_types_labels];

end