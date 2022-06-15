function f_dv_save_mat_data(app)
disp('Saving....');
fpath = app.matdatapathEditField.Value;

% remove things in case
% for n_dset = 1:numel(app.data.stats)
%     if isfield(app.data.stats{n_dset}.stat_params, 'cdata')
%         app.data.stats{n_dset}.stat_params = rmfield(app.data.stats{n_dset}.stat_params, 'cdata');
%     end
%     if isfield(app.data.stats{n_dset}.stat_params, 'ddata')
%         app.data.stats{n_dset}.stat_params = rmfield(app.data.stats{n_dset}.stat_params, 'ddata');
%     end
% end

var_list = [app.gui_ops.save_var_list_pl, app.gui_ops.save_var_list];

idx1 = strcmpi(app.data.Properties.VariableNames, 'dset_name_full');

for n_var = 1:numel(var_list)
    idx1 = idx1 + strcmpi(app.data.Properties.VariableNames, var_list{n_var});
end

idx1 = logical(idx1);

data_computed = app.data(:,idx1);

save(fpath, 'data_computed', '-v7.3');

disp('Done....');
end