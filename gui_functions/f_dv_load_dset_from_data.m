function f_dv_load_dset_from_data(app)

current_dset = app.DatasetDropDown.Value;

idx1 = strcmpi(app.data.experiment, current_dset);
app.ddata = app.data(idx1,:);

app.mplSpinner.Limits = [1, app.ddata.num_planes];

f_dv_update_dset_info(app);

end