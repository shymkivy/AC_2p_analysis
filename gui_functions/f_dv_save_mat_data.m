function f_dv_save_mat_data(app)
disp('Saving....');
fpath = app.matdatapathEditField.Value;

var_list = app.gui_ops.save_var_list;

idx1 = strcmpi(app.data.Properties.VariableNames, 'experiment');

for n_var = 1:numel(var_list)
    idx1 = idx1 + strcmpi(app.data.Properties.VariableNames, var_list{n_var});
end

idx1 = logical(idx1);

data_computed = app.data(:,idx1);

save(fpath, 'data_computed', '-v7.3');

disp('Done....');
end