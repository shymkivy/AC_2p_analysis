function f_dv_load_params(app)

app.BaselineEditField.Value = app.ops.trial_window(1);
app.RespEditField.Value = app.ops.trial_window(2);

app.regiontoplotDropDown.Items = [{'all'} app.ops.regions_to_analyze];

end