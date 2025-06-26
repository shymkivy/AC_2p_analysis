function f_dv_plot_reg_rois(app)

if strcmpi(app.regdatatouseDropDown.Value, 'caiman reg')
    source = 'register_roi_caiman_load';
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
    dset1 = data3(all_combs(n_comb,1),:);
    dset2 = data3(all_combs(n_comb,2),:);

    reg_data1 = dset1.(source){n_pl};
    reg_data2 = dset2.(source){n_pl};

    if ~isempty(reg_data1)
    
        reg_mat = [reg_data1.reg_cell_idx, reg_data2.reg_cell_idx];
        keep_cells = ~(sum(isnan(reg_mat),2));
        reg_mat2 = reg_mat(keep_cells,:);
        num_cells = size(reg_mat2,1);

        est1 = dset1.OA_data{n_pl}.est;
        est2 = dset2.OA_data{n_pl}.est;

        A1 = reshape(full(est1.A(:,reg_mat2(:,1))), [est1.dims(1), est1.dims(2), num_cells]);
        A2 = reshape(full(est2.A(:,reg_mat2(:,2))), [est2.dims(1), est2.dims(2), num_cells]);
        
        A1_im = sum(A1,3)';
        A2_im = sum(A2,3)';
        norm_fac = max(max([A1_im, A2_im]))/1.5;

        rgb_im = cat(3, A1_im/norm_fac, A2_im/norm_fac, A1_im/norm_fac);

        figure; imagesc(rgb_im); title(sprintf('%s im%d_%s(m) vs im%d_%s(g), pl%d; %d cells; %s', ddata.mouse_id{1}, dset1.im_num, dset1.dset_name{1}, dset2.im_num , dset2.dset_name{1}, n_pl, num_cells, app.regdatatouseDropDown.Value), 'interpreter', 'none')
    else
        fprintf('no data\n')
    end
end


end