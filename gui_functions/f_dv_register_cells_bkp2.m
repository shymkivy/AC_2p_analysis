function f_dv_register_cells_bkp2(app)

ddata = app.ddata;

%[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

idx1 = logical(strcmpi(ddata.mouse_id, app.data.mouse_id).*(ddata.FOV_num == app.data.FOV_num));

data2 = app.data(idx1,:);

input1 = app.regdsetuseEditField.Value;
if ~isempty(input1)

    dsets_use = zeros(numel(input1),1);
    for n_char = 1:numel(input1)
        dsets_use(n_char) = str2double(input1(n_char));
    end
    data2 = data2(dsets_use,:);
end

num_dsets = numel(data2.mouse_id);
num_planes = ddata.num_planes;
% 
% bkg_all = cell(num_dsets,1);
% 
% for n_dset = 1:num_dsets
%     est = data3(n_dset,:).OA_data{1}.est;
%     
%     bkg = reshape(mean(est.f)*est.b, [256 256]);
%     
%     A_orig = est.A * (mean(est.C,2) + mean(est.YrA,2));
%     A_orig = reshape(A_orig, est.dims');
%     
%     ave_im = bkg + A_orig;
%     
%     figure; 
%     imagesc(ave_im); axis equal tight
%     title(data3(n_dset,:).dset_name_full, 'interpreter', 'none');
%     
%     bkg_all{n_dset} = ave_im;
% end
% 
% [d1, d2] = size(bkg_all{1});
% 
% shifts_zy_all = cell(num_dsets,1);
% shifts_zy_all{1} = [0 0];
% targ = bkg_all{1};
% bkg_all_sf = bkg_all;
% %mc_params.zero_sides = 0;
% for n_dset = 2:num_dsets
%     shifts_zy_all{n_dset} = f_suite2p_reg_compute(bkg_all{n_dset}, targ);
%     bkg_all_sf{n_dset} = f_suite2p_reg_apply(bkg_all{n_dset}, shifts_zy_all{n_dset});
%     %bkg_all_sf{n_idx} = f_mc_apply_frame_shift(bkg_all{n_idx}, shifts_zy_all{n_idx}, mc_params);
% end
% 
% 
% all_pairs = nchoosek(1:num_dsets, 2);
% im3 = zeros(d1, d2, 3);
% for n_pair = 1:size(all_pairs,1)
%     im3(:,:,1) = bkg_all_sf{all_pairs(n_pair,1)}/max(bkg_all_sf{all_pairs(n_pair,1)}(:))*2;
%     im3(:,:,2) = bkg_all_sf{all_pairs(n_pair,2)}/max(bkg_all_sf{all_pairs(n_pair,2)}(:))*2;
%     figure; imagesc(im3); title(sprintf('pairs %d and %d', all_pairs(n_pair,1), all_pairs(n_pair,2)))
% end
% 
% imsh = f_mc_apply_frame_shift(bkg_all_sf{3}, [10 0], mc_params);
% 
% shifts_zy = f_mc_compute_frame_shift(targ, imsh);
% im_corr = f_mc_apply_frame_shift(imsh, shifts_zy, mc_params);
% 
% figure; imagesc(targ)
% figure; imagesc(imsh)
% figure; imagesc(im_corr)
%dsall= f_suite2p_reg_compute(bkg_all{2}, bkg_all{1});
%Y_reg = uint16(f_suite2p_reg_apply(Y_reg, dsall{n_iter}));

%%

A_all = cell(num_planes, num_dsets);
A_all3 = cell(num_planes, num_dsets);
num_cells = zeros(num_planes, num_dsets);
for n_dset = 1:num_dsets
    for n_pl = 1:num_planes
        data3 = data2(n_dset,:);
        est = data3.OA_data{n_pl}.est;
        proc = data3.OA_data{n_pl}.proc;
        num_cells(n_pl, n_dset) = sum(proc.comp_accepted);
        A_full = reshape(full(est.A(:,logical(proc.comp_accepted))), est.dims(1), est.dims(2), []);
        A_flat = sum(A_full,3);
        %figure; imagesc(A_flat);
        %title(num2str(fov_idx(n_idx)));
        A_all3{n_pl, n_dset} = A_full;
        A_all{n_pl, n_dset} = A_flat/max(A_flat(:))*2;
    end
end


throw_thresh = 0.5;
match_stats_id_all = cell(num_planes,num_dsets);
num_cells_common = zeros(num_planes,num_dsets);
match_stats_scores_all = cell(num_planes,num_dsets);
match_stats_max_scores_all = cell(num_planes,num_dsets);
dsets_all = 1:num_dsets;
for n_pl = 1:num_planes
    for n_dset = 1:num_dsets
        dsets_use = dsets_all;
        dsets_use(n_dset) = [];
        % first dset 1
        A_full = A_all3{n_pl, n_dset};
        num_cells1 = num_cells(n_pl, n_dset);
        match_stats_id = zeros(num_cells1, num_dsets);
        match_stats_num = zeros(num_cells1, num_dsets);
        match_stats_vals = cell(num_cells1, num_dsets);
        match_stats_max_val = zeros(num_cells1, num_dsets);
        match_stats_id(:,n_dset) = 1:num_cells1;
        for dset_idx = 1:(num_dsets-1)
            % compute scores of overlapping rois
            n_dset2 = dsets_use(dset_idx);
            A_full2 = A_all3{n_pl, n_dset2};
            for n_cell = 1:num_cells1
                score1 = squeeze(sum(sum(A_full(:,:,n_cell) .* A_full2,1),2));
                match_stats_num(n_cell, n_dset2) = sum(logical(score1));
                if match_stats_num(n_cell, n_dset2)
                    [match_stats_max_val(n_cell, n_dset2), match_stats_id(n_cell, n_dset2)] = max(score1);
                    match_stats_vals{n_cell, n_dset2} = score1(logical(score1));
                end
            end
            for n_cell = 1:num_cells1
                % for multiple overlaps pick highest
                idx1 = match_stats_id(n_cell,n_dset2) == match_stats_id(:,n_dset2);
                if sum(idx1) > 1
                    id1 = find(idx1);
                    [max_val1 , id2] = max(match_stats_max_val(id1, n_dset2));
                    throw = id1;
                    if max_val1 > throw_thresh
                        throw(id2) = [];
                    end
                    match_stats_id(throw,n_dset2) = 0;
                else
                    % discard low scores
                    if match_stats_max_val(n_cell, n_dset2) < throw_thresh
                        match_stats_id(n_cell,n_dset2) = 0;
                    end
                end
            end
        end

        idx_all = logical(prod(logical(match_stats_id),2));
        match_stats_id2 = match_stats_id(idx_all,:);
        
        match_stats_scores_all{n_pl, n_dset} = match_stats_vals;
        match_stats_max_scores_all{n_pl, n_dset} = match_stats_max_val;
        match_stats_id_all{n_pl, n_dset} = match_stats_id;
        num_cells_common(n_pl, n_dset) = size(match_stats_id2,1);
    end
end

n_pl = 1;
dset_idx = cell(num_dsets,1);
for n_dset = 1:num_dsets
    dset_idx{n_dset} = ones(num_cells(n_pl, n_dset),1)*n_dset;
end
dset_idx2 = cat(1, dset_idx{:});

% remove duplicates
temp_match = cat(1,match_stats_id_all{1,:});
temp_match2 = [(1:size(temp_match,1))', dset_idx2, temp_match];
num_rows = size(temp_match2,1);
n_row = 1;
while n_row <= num_rows
    score1 = sum(temp_match2(n_row,3:end) == temp_match2(:,3:end),2);
    score1(n_row) = 0;
    rem_idx = score1 == 4;
    temp_match2(rem_idx,:) = [];
    num_rows = size(temp_match2,1);
    n_row = n_row + 1;
end


A_rgb_all = cell(num_planes,1);
for n_pl = 1:num_planes
    temp_A1 = sum(A_all3{n_pl, 1}(:,:,match_stats_id_all{n_pl}(:,1)),3);
    temp_A2 = sum(A_all3{n_pl, 2}(:,:,match_stats_id_all{n_pl}(:,2)),3);
    temp_A3 = sum(A_all3{n_pl, 3}(:,:,match_stats_id_all{n_pl}(:,3)),3);
    A_rgb_all{n_pl} = cat(3,temp_A1 , temp_A2, temp_A3);
end

A_rgb_all2 = cat(2,A_rgb_all{:});

figure; imagesc(A_rgb_all2*5); axis equal tight
title(sprintf('%s; %d cells in common', ddata.dset_name_full{1}, sum(num_cells_common)), 'interpreter', 'none')

reg_out.match_stats_id = match_stats_id_all;
reg_out.num_cells_common = num_cells_common;

idx1 = logical(strcmpi(ddata.mouse_id, app.data.mouse_id).*(ddata.FOV_num == app.data.FOV_num));
idx2 = find(idx1, 1);
app.data(idx2,:).registration{1} = reg_out;


end