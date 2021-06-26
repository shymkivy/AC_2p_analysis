function f_dv_load_dset_from_data(app)

current_dset = app.DatasetDropDown.Value;
idx1 = strcmpi(app.data.experiment, current_dset);
app.current_data_idx = idx1;

%%
app.ddata = app.data(idx1,:);

%%
f_dv_load_ops(app);
f_dv_update_dset_info(app);
f_dv_update_A(app);
f_dv_update_cell(app);
f_dv_initialize_contours(app);
f_dv_set_contorus(app);

end