function f_dv_lick_tuning_similarity(app)


ddata = app.ddata;
idx1 = logical(strcmpi(ddata.mouse_id, app.data.mouse_id).*(ddata.FOV_num == app.data.FOV_num));
data2 = app.data(idx1,:);


stim_tone_idx = ddata.proc_data{1}.stim_params.ops.dev_tone_list;

stim_ops = data2(3,:).proc_data{1}.stim_params.ops;


stim_ops.MMN_patterns(stim_ops.paradigm_MMN_pattern(2),:)

num_planes = ddata.num_planes;


match_stats_id = data2.registration{1}.match_stats_id;

num_dsets = size(match_stats_id{1},2);

vec_all = cell(num_dsets,1);

n_dset = 1;
vec_all2 = cell(num_planes,1);
for n_pl = 1:num_planes
    stats2 = data2.stats_within{n_dset,n_pl};
    trial_data_mean_z = stats2.stim_trial_data_mean_z;
    peak_vals = max(trial_data_mean_z, [], 2);
    vec_all2{n_pl} = peak_vals(match_stats_id{n_pl}(:,n_dset));
end    
vec_all{1} = cat(1, vec_all2{:});


n_dset = 3;
vec_all2 = cell(num_planes,1);
for n_pl = 1:num_planes
    stats1 = data2.stats{n_dset,n_pl};
    match_stats_id2 = match_stats_id{n_pl};
    peak_vals = stats1.peak_val_all(match_stats_id2(:,n_dset), stim_tone_idx);
    
    vec_all2{n_pl} = peak_vals;
end    
vec_all{2} = cat(1, vec_all2{:});


n_dset = 3;
vec_all2 = cell(num_planes,1);
for n_pl = 1:num_planes
    stats1 = data2.stats{n_dset,n_pl};
    match_stats_id2 = match_stats_id{n_pl};
    peak_vals = stats1.peak_val_all(match_stats_id2(:,n_dset), 20);
    
    vec_all2{n_pl} = peak_vals;
end    
vec_all{3} = cat(1, vec_all2{:});


vec_all3 = reshape(cat(3, vec_all{:}), [], num_dsets);

vec_all3 = vec_all3./std(vec_all3);


si_cos = 1 - pdist2(vec_all3', vec_all3', 'cosine');

si_corr = 1 - pdist2(vec_all3', vec_all3', 'correlation');

figure; imagesc(si_cos); axis equal tight;
title(sprintf('%s; population trial ave cosine similarity', ddata.mouse_id{1}), 'interpreter', 'none');

figure; imagesc(si_corr); axis equal tight;
title(sprintf('%s; population trial ave correlation', ddata.mouse_id{1}), 'interpreter', 'none');


% 
% tn_all_sel = f_dv_get_trial_number(app);
% %tt_all = app.ops.context_types_all(tn_all_sel)';
% 
% cmap_all = app.ops.context_types_all_colors;
% 
% cmap_all(18,:,:) = [1 1 0];
% cmap_all(28,:,:) = [1 1 0];

%cmap1 = jet(10);
%cmap2 = reshape(cmap1, 10, 1, 3);

cmap_gray = [.6 .6 .6];
cmap_gray2 = reshape(cmap_gray, 1, 1, 3);


pl_im_all = cell(ddata.num_planes,1);
for n_pl = 1:ddata.num_planes
    dims = ddata.OA_data{n_pl}.est.dims;
    comp_accepted = ddata.OA_data{n_pl}.proc.comp_accepted;
    num_cells = sum(comp_accepted);
    stats1 = ddata.stats{n_pl};
    
    A = reshape(full(ddata.OA_data{n_pl}.est.A(:,comp_accepted)), dims(1), dims(2), num_cells);
    A_n = A/max(A(:));

    A_col = zeros(dims(1), dims(2), 3);

    for n_cell = 1:num_cells
        temp_A = A_n(:,:,n_cell);
        %temp_A = temp_A/max(temp_A(:));
        temp_A2 = repmat(temp_A, 1, 1, 3);
        if sum(stats1.resp_cells_peak(n_cell,1:10))
            [~, n_freq] = max(stats1.peak_val_all(n_cell, tn_all_sel));
            temp_A_col = temp_A2.*cmap_all(tn_all_sel(n_freq),:,:);
        else
            temp_A_col = temp_A2.*cmap_gray2;
        end
        A_col = A_col + temp_A_col;
    end

    pl_im_all{n_pl} = A_col*2;
end
pl_im_all_flat = cat(2, pl_im_all{:});

figure; imagesc(pl_im_all_flat); axis equal tight;
title(sprintf('%s; trias [%s]', ddata.dset_name_full{1}, num2str(tn_all_sel)), 'interpreter', 'none');

figure; imagesc(cmap_all(tn_all_sel,:,:))


end