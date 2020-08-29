function [data, ops] = f_mpl_load_data(ops)

disp('Loading data...');

% gather the file names of specified data
dset_index = cell(numel(ops.regions_to_analyze),1);
dset_total_count = 0;
for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    temp_file_names = ops.file_names.(cond_name);
    data.(cond_name).num_dsets = numel(temp_file_names);
    data.(cond_name).num_planes = zeros(data.(cond_name).num_dsets,1);
    for n_dset = 1:numel(temp_file_names)
        temp_proc = dir([ops.file_dir, ['\' cond_name '*' ops.paradigm_type '*' temp_file_names{n_dset} '*' ops.processed_data_tag '.mat']]);
        temp_OA = dir([ops.file_dir, ['\' cond_name '*' ops.paradigm_type '*' temp_file_names{n_dset} '*' ops.OA_output_tag '.mat']]);
        
        if isempty(temp_proc)
            error(['S12 processed data file missing in ' cond_name ' ' temp_file_names{n_dset}])
        end
        if isempty(temp_OA)
            error(['Sorted data file missing in ' cond_name ' ' temp_file_names{n_dset}])
        end
        
        for n_pl = 1:numel(temp_proc)
            data.(cond_name).file_names_proc_data{n_dset,n_pl} = temp_proc(n_pl).name;
        end
        for n_pl = 1:numel(temp_OA)
            data.(cond_name).file_names_OA_results{n_dset,n_pl} = temp_OA(n_pl).name;
        end
        data.(cond_name).num_planes(n_dset,1) = numel(temp_OA);
    end
    dset_index{n_cond} = (1:data.(cond_name).num_dsets)'+dset_total_count;
    dset_total_count = dset_index{n_cond}(end);
end
ops.dset_index = dset_index;
ops.dset_total_count = dset_total_count;

%% load
if ops.waitbar
    wb = f_waitbar_initialize([], 'Loading data...');
end
for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    for n_dset = 1:data.(cond_name).num_dsets
        
        % load proc data
        temp_load = load([ops.file_dir '\' data.(cond_name).file_names_proc_data{n_dset}]);
        data.(cond_name).proc_data{n_dset,1} = temp_load.data;
        data.(cond_name).proc_ops{n_dset,1} = temp_load.ops;
        
        % load OA data
        for n_pl = 1:size(data.(cond_name).file_names_OA_results,2)
            data.(cond_name).OA_data{n_dset,n_pl} = load([ops.file_dir '\' data.(cond_name).file_names_OA_results{n_dset,n_pl}]);
        end
        
        if ops.waitbar
            f_waitbar_update(wb, dset_index{n_cond}(n_dset)/dset_total_count, sprintf('Loading %s, dset %d', cond_name, n_dset));
        end
    end
end
if ops.waitbar
    f_waitbar_close(wb);
end
end