function f_dv_load_data_mat(app)
% unuseed

fpath = app.s31datapathEditField.Value;

disp('Loading data...')
data_load = load(fpath);

app.data = data_load.data;
app.ops = data_load.ops;

%%
f_dv_initialize_post_load(app);

end