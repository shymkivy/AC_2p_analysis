function f_dv_load_data_button(app)
%% load preprocessing parameters
ops = app.ops;
ops.app = app;

ops.experiment_type = app.ExperimentDropDown.Value;
ops.load_mat_data = app.automatCheckBox.Value;
ops.load_reg_data = app.autoregCheckBox.Value;
ops.params = f_dv_gather_params(app);

%% List of files to load
[app.data, ops, reg_struct] = f_load_data(ops);

if isfield(reg_struct, 'reg_data')
    app.reg_data = reg_struct.reg_data;
    app.border_coords = reg_struct.border_coords;
end

%%
ops = rmfield(ops, 'app');
app.ops = ops;

%%
f_dv_initialize_post_load(app);
f_dv_load_dset_from_data(app);
end