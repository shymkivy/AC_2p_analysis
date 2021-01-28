function f_fov_registration(data, ops)
%% for each dset plot things

reg_data_path = 'C:\Users\ys2605\Desktop\stuff\register_2p_to_wf\reg_save.mat';
reg_data = load(reg_data_path);

mouse_names = unique(ops.AC_data.mouse_tag);

areas_all = fieldnames(data);

for n_mous = 1:numel(mouse_names)
    dset_idx = strcmpi(mouse_names{n_mous},ops.AC_data.mouse_tag);
    temp_AC_data = ops.AC_data(dset_idx,:);
   
    temp_field = areas_all(strcmpi(areas_all, temp_AC_data.area(1)));
    
    cdata = data.(temp_field{1});
    
    temp_AC_data.experiment(1)
    
end

for n_cond = 1:numel(ops.regions_to_analyze)
    cond_name = ops.regions_to_analyze{n_cond};
    cdata = data.(cond_name);
    
    for n_dset = 1:cdata.num_dsets
    end
end


end