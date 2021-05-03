function [data, ops] = f_mpl_load_data(ops)

disp('Loading data...');
data = ops.AC_data;

% gather the file names of specified data
num_dset = numel(data.area);
if ops.waitbar
    wb = f_waitbar_initialize([], 'Loading data...');
end
for n_dset = 1:num_dset
    fname = data.experiment{n_dset};
    
    % load proc data
    temp_proc = dir([ops.file_dir, ['\' '*' fname '*' ops.processed_data_tag '.mat']]);
    if isempty(temp_proc)
        error(['S12 processed data file missing: ' fname])
    end
    temp_load = load([ops.file_dir '\' temp_proc.name]);
    data.proc_data{n_dset} = temp_load.data;
    data.proc_ops{n_dset} = temp_load.ops;

    % load OA
    temp_OA = dir([ops.file_dir, ['\' '*' fname '*' ops.OA_output_tag '.mat']]);
    if isempty(temp_OA)
        error(['S12 sorted OA data file missing: ' fname])
    end
    data.num_planes(n_dset) = numel(temp_OA);
    for n_pl = 1:data.num_planes(n_dset)
        data.OA_data{n_dset, n_pl} = load([ops.file_dir '\' temp_OA(n_pl).name]);
    end
    
    if ops.waitbar
        f_waitbar_update(wb, n_dset/num_dset, sprintf('Loading %d/%d',n_dset,num_dset));
    end
end

if ops.waitbar
    f_waitbar_close(wb);
end
    
end