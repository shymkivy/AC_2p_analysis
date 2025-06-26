function data = f_dv_load_mat_data(data, ops)

disp('Loading mat data');

fpath = ops.fpath_mat_data;

if exist(fpath, 'file') == 2
    data1 = load(fpath);
    data1 = data1.data_computed;
    
    num_dsets = size(data,1);
    
    var_list = [ops.save_var_list_pl, ops.save_var_list];
    
    for n_dset = 1:num_dsets
        idx1 = strcmpi(data(n_dset,:).dset_name_full, data1.dset_name_full);
        if sum(idx1)
            for n_var = 1:numel(var_list)
                var1 = var_list{n_var};
                if sum(strcmpi(data1.Properties.VariableNames, var1))
                    for n_pl = 1:numel(data1.(var1)(idx1,:))
                        if ~isempty(data1(idx1,:).(var1){n_pl})
                            data(n_dset,:).(var1){n_pl} = data1(idx1,:).(var1){n_pl};
                        end
                    end
                end
            end
        end
    end
    
    disp('Done');
else
    disp('Mat data file not available')
end

end