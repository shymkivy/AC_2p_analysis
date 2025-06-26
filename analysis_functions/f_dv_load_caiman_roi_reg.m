function f_dv_load_caiman_roi_reg(app)
disp('Loading caiman reg data...');
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
                mpl_tag = sprintf('_mpl%d_pl%d', num_planes, n_pl);
            else
                mpl_tag = '';
            end
            
            file1 = dir(sprintf('%s\\*%s*fov%d*%s_registration_cmnf.mat',data_dir, mouse_all{n_ms}, fov1, mpl_tag));
            
            if numel(file1)
                reg_data_load = load([data_dir '\' file1.name]);
                
                reg_mat = reg_data_load.reg_out{2} + 1;
                num_dsets = numel(data3.dset_name_full);
                 
                % first for every dset find reg data
                dset_reg = cell(num_dsets,1);
                for n_dset2 = 1:num_dsets
                    fname_dset = [data3.dset_name_full{n_dset2} mpl_tag];
                    for n_fl = 1:size(reg_mat,2)
                        fname2 = strtrim(reg_data_load.fname_list(n_fl,:));
                        if strcmpi(fname_dset, fname2)
                            dset_reg{n_dset2} = reg_mat(:,n_fl);
                        end
                    end
                end
                dset_reg2 = cat(2,dset_reg{:});
                keep_cell = ~(sum(isnan(dset_reg2),2) == size(dset_reg2,2));
                
                reg_data.fov_cell_idx = (1:sum(keep_cell))';
                
                for n_dset2 = 1:num_dsets
                    reg_data2 = reg_data;
                    reg_data2.reg_cell_idx = dset_reg{n_dset2}(keep_cell);
                    idx1 = data.idx == data3.idx(n_dset2);
                    app.data.register_roi_caiman_load{idx1, n_pl} = reg_data2;
                end
            end
        end
    end
end
disp('Done');

end