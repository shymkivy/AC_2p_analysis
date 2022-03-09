function [data, ops] = f_mpl_load_data(ops)

disp('Loading data...');
data = ops.AC_data;

% gather the file names of specified data
num_dset = numel(data.area);
if ops.waitbar
    wb = f_waitbar_initialize([], 'Loading data...');
end

has_data = false(num_dset,1);
for n_dset = 1:num_dset
    ddata = data(n_dset,:);   
    fname_tag = ddata.dset_name{1};
    date_tag = ddata.mouse_tag{1};
    
    % load proc data
    temp_proc = dir([ops.file_dir, ['\' '*' fname_tag '*' date_tag '*' ops.processed_data_tag '.mat']]);
    if ~isempty(temp_proc)
        temp_load = load([ops.file_dir '\' temp_proc.name]);
        data.proc_data{n_dset} = temp_load.data;
        data.proc_ops{n_dset} = temp_load.ops;
    end
    

    % load OA
    temp_OA = dir([ops.file_dir, ['\' '*' fname_tag '*' date_tag '*' ops.OA_output_tag '.mat']]);
    if ~isempty(temp_OA)
        has_data(n_dset) = 1;
        data.num_planes(n_dset) = numel(temp_OA);
        for n_pl = 1:data.num_planes(n_dset)
            data.OA_data{n_dset, n_pl} = load([ops.file_dir '\' temp_OA(n_pl).name]);
        end
    end
    
    
    if ops.waitbar
        f_waitbar_update(wb, n_dset/num_dset, sprintf('Loading %d/%d',n_dset,num_dset));
    end
end

if ops.waitbar
    f_waitbar_close(wb);
end
data = data(has_data,:);

end