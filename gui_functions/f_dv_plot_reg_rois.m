function f_dv_plot_reg_rois(app)

if strcmpi(app.regdatatouseDropDown.Value, 'caiman reg')
    source = 'registration_caiman';
elseif strcmpi(app.regdatatouseDropDown.Value, 'gui reg')
    source = 'register_roi';
end

n_pl = app.mplSpinner.Value;

data = app.data;
ddata = data(app.current_data_idx,:);
data2 = data(strcmpi(data.mouse_id, ddata.mouse_id),:);
data3 = data2(data2.FOV_num == ddata.FOV_num,:);

num_dsets = numel(data3.mouse_id);

all_combs = nchoosek(1:num_dsets,2);

for n_comb = 1:size(all_combs,1)
    dset1 = all_combs(n_comb,1);
    dset2 = all_combs(n_comb,2);
    reg_data1 = data3.(source){dset1, n_pl};
    reg_data2 = data3.(source){dset2, n_pl};
    
    if ~isempty(reg_data1)
    
        est1 = data3.OA_data{dset1, n_pl}.est;
        est2 = data3.OA_data{dset2, n_pl}.est;

        proc1 = data3.OA_data{dset1, n_pl}.proc;
        proc2 = data3.OA_data{dset2, n_pl}.proc;

        reg_mat = [reg_data1.reg_cell_idx, reg_data2.reg_cell_idx];
        keep_cells = ~(sum(isnan(reg_mat),2));
        reg_mat2 = reg_mat(keep_cells,:);


        if strcmpi(app.regdatatouseDropDown.Value, 'caiman reg')
            
            cell_idx_cut1 = (1:sum(proc1.idx_cell_acc))';
            idx_cell_lab1 = double(proc1.idx_cell_acc);
            idx_cell_lab1(proc1.idx_cell_acc) = cell_idx_cut1;

            cell_idx_cut2 = (1:sum(proc2.idx_cell_acc))';
            idx_cell_lab2 = double(proc2.idx_cell_acc);
            idx_cell_lab2(proc2.idx_cell_acc) = cell_idx_cut2;
            
            reg_mat3 = [idx_cell_lab1(reg_mat2(:,1)), idx_cell_lab2(reg_mat2(:,2))];
            keep_cells2 = logical(prod(logical(reg_mat3),2));
            reg_mat4 = reg_mat3(keep_cells2,:);
        elseif strcmpi(app.regdatatouseDropDown.Value, 'gui reg')
            reg_mat4 = reg_mat2;
        end

        num_cells = size(reg_mat4,1);

        A1 = reshape(full(est1.A(:,reg_mat4(:,1))), [est1.dims(1), est1.dims(2), num_cells]);
        A2 = reshape(full(est2.A(:,reg_mat4(:,2))), [est2.dims(1), est2.dims(2), num_cells]);
        
        cells1 = find(proc1.idx_cell_acc);
        cells2 = find(proc2.idx_cell_acc);

        figure;
        imagesc(reshape(full(est1.A(:,cells1(3))), [est2.dims(1), est2.dims(2), 1]))
        
        figure;
        imagesc(reshape(full(est2.A(:,cells2(3))), [est2.dims(1), est2.dims(2), 1]))

        A1_im = sum(A1,3)';
        A2_im = sum(A2,3)';
        norm_fac = max(max([A1_im, A2_im]))/1.5;

        rgb_im = cat(3, A1_im/norm_fac, A2_im/norm_fac, A1_im/norm_fac);

        figure; imagesc(rgb_im); title(sprintf('%s im%d_%s(m) vs im%d_%s(g), pl%d; %s', ddata.mouse_id{1}, data3.im_num(dset1), data3.dset_name{dset1}, data3.im_num(dset2) , data3.dset_name{dset2}, n_pl, app.regdatatouseDropDown.Value), 'interpreter', 'none')
    else
        fprintf('no data\n')
    end
end


end