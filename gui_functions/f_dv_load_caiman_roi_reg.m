function f_dv_load_caiman_roi_reg(app)

data = app.data;
data_dir = app.ops.file_dir;

mouse_all = unique(data.mouse_id, 'stable');

for n_ms = 1:numel(mouse_all)
    data2 = data(strcmpi(data.mouse_id, mouse_all{n_ms}),:);
    fov_all = unique(data2.FOV_num, 'stable');
    for n_fov = 1:numel(fov_all)
        fov1 = fov_all(n_fov);
        data3 = data2(data2.FOV_num == fov1,:);
        
        num_planes = data3.num_planes(1);
        
        for n_pl = 1:num_planes
            if num_planes > 1
                file1 = dir(sprintf('%s\\*%s*fov%d*pl%d*registration_cmnf.mat',data_dir, mouse_all{n_ms}, fov1, n_pl));
            else
                file1 = dir(sprintf('%s\\*%s*fov%d*registration_cmnf.mat',data_dir, mouse_all{n_ms}, fov1));
            end
            
            if numel(file1)
                reg_data = load([data_dir '\' file1.name]);
                idx1 = data.idx == data3.idx(1);
                app.data.registration_caiman{idx1, n_pl} = reg_data;
            end
        end
    end
end

end