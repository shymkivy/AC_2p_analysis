function [data, ops] = f_mpl_load_data2(ops)

disp('Loading data...');

% gather the file names of specified data
for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    temp_file_names = ops.file_names.(cond_name);
    %data.(cond_name).num_dsets = numel(temp_file_names);
    %data.(cond_name).num_planes = zeros(data.(cond_name).num_dsets,1);
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
            data.(cond_name)(n_dset).file_names_proc_data{n_pl,1} = temp_proc(n_pl).name;
        end
        for n_pl = 1:numel(temp_OA)
            data.(cond_name)(n_dset).file_names_OA_results{n_pl,1} = temp_OA(n_pl).name;
        end
        data.(cond_name)(n_dset).num_planes = numel(temp_OA);
    end
end

%ops.file_dir '\' 

wb = f_waitbar_initialize([], 'Loading data...');

for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    for n_dset = 1:numel(data.(cond_name))
        
        % load proc data
        temp_load = load([ops.file_dir '\' data.(cond_name)(n_dset).file_names_proc_data{:}]);
        data.(cond_name)(n_dset).proc_data = temp_load.data;
        data.(cond_name)(n_dset).proc_ops = temp_load.ops;
        
        % load OA data
        for n_pl = 1:numel(data.(cond_name)(n_dset).file_names_OA_results)
            data.(cond_name)(n_dset).OA_data{n_pl,1} = load([ops.file_dir '\' data.(cond_name)(n_dset).file_names_OA_results{n_pl}]);
        end
        
        f_waitbar_update(wb, (n_dset/numel(data.(cond_name)) + (n_cond-1))/numel(ops.regions_to_analyze));
        
    end
end

f_waitbar_close(wb);

end