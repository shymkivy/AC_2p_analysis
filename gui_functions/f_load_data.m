function [data, ops, reg_struct] = f_load_data(ops)

ops = f_AC_files_import(ops);

ops = f_process_ops(ops);

[data, ops] = f_load_datasets(ops);

data = f_dv_preprocess_data(data, ops);

data = f_dv_compute_cdata_all(data, ops);

% load more
if ops.load_mat_data
    data = f_dv_load_mat_data(data, ops);
end

if ops.load_reg_data
    [reg_struct, data] = f_dv_load_reg_data(data, ops);
end

end