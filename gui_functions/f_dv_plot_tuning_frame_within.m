function f_dv_plot_tuning_frame_within(app)


ddata = app.ddata;

%idx1 = logical(strcmpi(ddata.mouse_id, app.data.mouse_id).*(ddata.FOV_num == app.data.FOV_num));
%data2 = app.data(idx1,:);
stim_tone_idx = ddata.proc_data{1,1}.stim_params.ops.dev_tone_list;

% 
% tn_all_sel = f_dv_get_trial_number(app);
% %tt_all = app.ops.context_types_all(tn_all_sel)';
% 
cmap_all = app.ops.context_types_all_colors; 
cmap_all(18,:,:) = [1 1 0];
cmap_all(28,:,:) = [1 1 0];

%cmap1 = jet(10);
%cmap2 = reshape(cmap1, 10, 1, 3);

cmap_gray = [.6 .6 .6];
cmap_gray2 = reshape(cmap_gray, 1, 1, 3);

z_thresh = app.TFwithinzthreshEditField.Value;

pl_im_all = cell(ddata.num_planes,1);
for n_pl = 1:ddata.num_planes
    dims = ddata.OA_data{n_pl}.est.dims;
    comp_accepted = ddata.OA_data{n_pl}.proc.comp_accepted;
    num_cells = sum(comp_accepted);
    stats1 = ddata.stats_within{n_pl};
    
    peak_resp = max(stats1.stim_trial_data_mean_z, [], 2);
    
    A = reshape(full(ddata.OA_data{n_pl}.est.A(:,comp_accepted)), dims(1), dims(2), num_cells);
    A_n = A/max(A(:));

    A_col = zeros(dims(1), dims(2), 3);

    for n_cell = 1:num_cells
        temp_A = A_n(:,:,n_cell);
        %temp_A = temp_A/max(temp_A(:));
        temp_A2 = repmat(temp_A, 1, 1, 3);
        if peak_resp(n_cell) > z_thresh
            temp_A_col = temp_A2.*cmap_all(stim_tone_idx,:,:)*peak_resp(n_cell)/max(peak_resp)*2;
        else
            temp_A_col = temp_A2.*cmap_gray2;
        end
        A_col = A_col + temp_A_col;
    end

    pl_im_all{n_pl} = A_col*2;
end
pl_im_all_flat = cat(2, pl_im_all{:});

figure; imagesc(pl_im_all_flat); axis equal tight;
title(sprintf('%s; tone resp cells, trial %d, z=%d', ddata.dset_name_full{1}, stim_tone_idx, z_thresh), 'interpreter', 'none');

figure; imagesc(cmap_all(stim_tone_idx,:,:))

end