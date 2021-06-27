function f_dv_load_mat_data(app)

fpath = app.matdatapathEditField.Value;

data1 = load(fpath);
data1 = data1.data_computed;

num_dsets = size(data1,1);

var_list = app.gui_ops.save_var_list;

for n_dset = 1:num_dsets
    idx1 = strcmpi(app.data.experiment, data1(n_dset,:).experiment);
    if sum(idx1)
        for n_var = 1:numel(var_list)
            var1 = var_list{n_var};
            if sum(strcmpi(data1.Properties.VariableNames, var1))
                for n_pl = 1:numel(data1.(var1)(n_dset))
                    if ~isempty(data1(n_dset,:).(var1){n_pl})
                        app.data(idx1,:).(var1){n_pl} = data1(n_dset,:).(var1){n_pl};
                    end
                end
            end
        end
    end
end

f_dv_initialize_post_load(app);

end