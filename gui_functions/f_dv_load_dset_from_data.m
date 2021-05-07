function f_dv_load_dset_from_data(app)

current_dset = app.DatasetDropDown.Value;

idx1 = strcmpi(app.data.experiment, current_dset);
app.ddata = app.data(idx1,:);

app.trialtypeDropDown.Items = [{'all'}; app.ops.context_types_labels];

f_dv_load_ops(app);
f_dv_update_dset_info(app);

end