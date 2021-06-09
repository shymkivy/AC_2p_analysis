function f_dv_load_data_mat(app)
% unuseed


fpath = app.s31datapathEditField.Value;

disp('Loading data...')
data_load = load(fpath);

app.data = data_load.data;
app.ops = data_load.ops;

app.DatasetDropDown.Items = app.data.experiment;

f_dv_load_dset_from_data(app);
disp('Loaded')
end