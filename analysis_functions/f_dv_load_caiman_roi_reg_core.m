function reg_out = f_dv_load_caiman_roi_reg_core(data, data_dir)

num_dsets = numel(data.mouse_id);
num_planes = data.num_planes(1);

mouse_id = data.mouse_id{1};
fov1 = data.FOV_num(1);

reg_out = cell(num_dsets, num_planes);
for n_pl = 1:num_planes
    if num_planes > 1
        mpl_tag = sprintf('_mpl%d_pl%d', num_planes, n_pl);
    else
        mpl_tag = '';
    end

    file1 = dir(sprintf('%s\\*%s*fov%d*%s_registration_cmnf.mat',data_dir, mouse_id, fov1, mpl_tag));

    if numel(file1)
        reg_data_load = load([data_dir '\' file1.name]);

        reg_mat = reg_data_load.reg_out{2} + 1;
        num_dsets = numel(data.dset_name_full);

        % first for every dset find reg data
        dset_reg = cell(num_dsets,1);
        for n_dset = 1:num_dsets
            fname_dset = [data.dset_name_full{n_dset} mpl_tag];
            for n_fl = 1:size(reg_mat,2)
                fname2 = strtrim(reg_data_load.fname_list(n_fl,:));
                if strcmpi(fname_dset, fname2)
                    dset_reg{n_dset} = reg_mat(:,n_fl);
                end
            end
        end
        dset_reg2 = cat(2,dset_reg{:});
        keep_cell = ~(sum(isnan(dset_reg2),2) == size(dset_reg2,2));

        reg_data.fov_cell_idx = (1:sum(keep_cell))';

        for n_dset = 1:num_dsets
            reg_data2 = reg_data;
            reg_data2.reg_cell_idx = dset_reg{n_dset}(keep_cell);
            reg_out{n_dset, n_pl} = reg_data2;
        end
    end
end


end