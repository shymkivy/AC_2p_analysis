function [data_all3, tt_all3, reg_id3, group_id3] = f_dv_decoder_gather_data(data0, params, ops)

tn_all = f_dv_get_trial_number(params);
[data, ~] = f_dv_get_data_by_mouse_selection(data0, params);
[cdata0, ~] = f_dv_get_new_cdata_stats(data(1,:), params);
[region_num, ~, ~] = f_dv_get_region_sel_val(params, ops);

[num_gr, num_tn] = size(tn_all);
num_dsets = size(data,1);
num_regions = size(region_num,1);

trial_window = [-1, 3];
[~, trial_frames] = f_dv_compute_window_t(trial_window, cdata0(1).volume_period);

reg_id = ones(num_dsets, num_gr, num_regions).*reshape((1:num_regions), 1, 1, num_regions);
group_id = ones(num_dsets, num_gr, num_regions).*reshape((1:num_gr), 1, num_gr, 1);
data_all = cell(num_dsets, num_gr, num_regions);
tt_all = cell(num_dsets, num_gr, num_regions);

for n_dset = 1:num_dsets
    % extracting firing rates data from data structure
    ddata = data(n_dset,:);
    [cdata, stats1] = f_dv_get_new_cdata_stats(ddata, params);
    firing_rate = cat(1,cdata.S_sm);
    
    mmn_freq = ddata.MMN_freq{1};
    stim_times = ddata.stim_frame_index{1};

    trial_types = ddata.trial_types{1};
    if ~isempty(mmn_freq)
        trial_types_ctx2 = f_dv_mark_tt_ctx(trial_types, mmn_freq, ops);
        trial_types_all = [trial_types, trial_types_ctx2];
    else
        trial_types_all = trial_types;
    end

    reg_cell_labels = f_dv_get_area_label(ddata, params, ops);

    for n_gr = 1:num_gr
        tn1 = tn_all(n_gr, :);
        resp_cells = f_dv_get_resp_vals_cells(stats1, tn1, params);
        resp_cells2 = logical(sum(resp_cells,2));
        for n_reg = 1:num_regions
            n_reg2 = region_num(n_reg,:);

            reg_cell_idx = logical(sum(reg_cell_labels == n_reg2,2));
            resp_reg_cell = and(resp_cells2, reg_cell_idx);
            
            num_cells2 = sum(resp_reg_cell);
            
            if num_cells2 > 10

                firing_rate2 = firing_rate(resp_reg_cell,:);
    
                trial_data_sort = f_get_stim_trig_resp(firing_rate2, stim_times, trial_frames);
                
                tt1 = reshape(ops.context_types_all(tn1), [1, 1, numel(tn1)]);
                idx1 = logical(sum(trial_types_all == tt1,3));
                num_tt2 = size(trial_types_all,2);
                trial_types2c = cell(num_tt2,1);
                trial_data_sort2c = cell(num_tt2,1);
                for n_r = 1:num_tt2
                    trial_types2c{n_r} = trial_types_all(idx1(:,n_r),n_r);
                    trial_data_sort2c{n_r} = trial_data_sort(:,:,idx1(:,n_r));
                end
    
                trial_types2 = cat(1, trial_types2c{:});
                trial_data_sort2 = cat(3,trial_data_sort2c{:});
                
                trial_counts = sum(trial_types2 == ops.context_types_all(tn1)',1);
                
                if strcmpi(params.trial_num_selection, 'all')
                    max_tr = max(trial_counts);
                elseif strcmpi(params.trial_num_selection, 'mean')
                    max_tr = round(mean(trial_counts));
                elseif strcmpi(params.trial_num_selection, 'median')
                    max_tr = round(defian(trial_counts));
                elseif strcmpi(params.trial_num_selection, 'min')
                    max_tr = round(min(trial_counts));
                end

                trial_types3c = cell(num_tn,1);
                trial_data_sort3c = cell(num_tn,1);
                for n_tt = 1:num_tn
                    tt2 = ops.context_types_all(tn1(n_tt));
                    tt_idx1 = trial_types2 == tt2;
                    if sum(tt_idx1) > max_tr
                        samp_tr1 = randsample(find(tt_idx1), max_tr);
                        trial_types3c{n_tt} = trial_types2(samp_tr1);
                        trial_data_sort3c{n_tt} = trial_data_sort2(:,:,samp_tr1);
                    else
                        trial_types3c{n_tt} = trial_types2(tt_idx1);
                        trial_data_sort3c{n_tt} = trial_data_sort2(:,:,tt_idx1);
                    end
                end

                trial_types3 = cat(1, trial_types3c{:});
                trial_data_sort3 = cat(3,trial_data_sort3c{:});
                
                data_all{n_dset, n_gr, n_reg} = trial_data_sort3;
                tt_all{n_dset, n_gr, n_reg} = trial_types3;
            end
        end
    end
end

data_all2 = reshape(data_all,[],1);
tt_all2 = reshape(tt_all,[],1);
reg_id2 = reshape(reg_id,[],1);
group_id2 = reshape(group_id,[],1);

num_dsets2 = numel(data_all2);
has_data = false(num_dsets2,1);
for n_data = 1:num_dsets2
    if ~isempty(data_all{n_data})
        has_data(n_data) = 1;
    end
end

data_all3 = data_all2(has_data);
tt_all3 = tt_all2(has_data);
reg_id3 = reg_id2(has_data);
group_id3 = group_id2(has_data);

end