function [features, sel_cells, area_labels, features_proc] = f_dv_get_feature(feature, tn_all, data, params, ops)

[data, ~] = f_dv_get_data_by_mouse_selection(data, params);

[region_num, ~, ~] = f_dv_get_region_sel_val(params, ops);
num_reg = size(region_num,1);

params2 = params;

num_dsets = size(data,1);
[num_gr, num_tn] = size(tn_all);

features = cell(num_dsets, num_gr, num_tn);
sel_cells = cell(num_dsets, num_gr, num_tn);
area_labels = cell(num_dsets, num_gr, 1);
features_proc = cell(num_dsets, num_gr, num_reg, num_tn);
for n_dset = 1:num_dsets
    stats1 = cat(1,data(n_dset,:).stats{params.planes});
    % first get responsive cells, then extract features
    
    for n_gr = 1:num_gr
        
        tn_all2 = tn_all(n_gr,:);
        params2.responsive_cells_type = 'peaks';
        [sel_cells_peak, ~, vals_peak] = f_dv_get_resp_vals_cells(stats1, tn_all2, params);
        params2.responsive_cells_type = 'onset';
        [sel_cells_onset, ~, vals_onset] = f_dv_get_resp_vals_cells(stats1, tn_all2, params);
        params2.responsive_cells_type = 'offset';
        [sel_cells_offset, ~, vals_offset] = f_dv_get_resp_vals_cells(stats1, tn_all2, params);
        
        % get area labels
        area_labels0 = f_dv_get_area_label(data(n_dset,:), params, ops);
        area_labels{n_dset, n_gr, 1} = area_labels0;
        for n_tn = 1:num_tn
            tn1 = tn_all2(n_tn);
            
            vals_peak2 = vals_peak(:, n_tn);
            vals_onset2 = vals_onset(:, n_tn);
            vals_offset2 = vals_offset(:, n_tn);

            resp_cells_peak1 = sel_cells_peak(:, n_tn);
            resp_cells_onset1 = sel_cells_onset(:, n_tn);
            resp_cells_offset1 = sel_cells_offset(:, n_tn);
            
            if strcmpi(feature, 'peak resp mag')
                sel_cells2 = resp_cells_peak1;
                features2 = vals_peak2;
            elseif strcmpi(feature, 'peak resp mag z')
                sel_cells2 = resp_cells_peak1;
                vals1 = vals_peak2;
                features2 = (vals1 - stats1.stat_trials_mean_mean)./stats1.stat_trials_mean_sem;
            elseif strcmpi(feature, 'peak loc')
                sel_cells2 = resp_cells_peak1;
                vals1 = cat(1, stats1.peak_loc);
                features2 = vals1(:,tn1);
            elseif strcmpi(feature, 'peak resp thresh')
                sel_cells2 = resp_cells_peak1;
                features2 = cat(1, stats1.peak_resp_thresh);
            elseif strcmpi(feature, 'stat trials mean sem')
                sel_cells2 = resp_cells_peak1;
                features2 = cat(1, stats1.stat_trials_mean_sem);
            elseif strcmpi(feature, 'onset mag')
                sel_cells2 = resp_cells_onset1;
                features2 = vals_onset2;
            elseif strcmpi(feature, 'onset mag z')
                sel_cells2 = resp_cells_onset1;
                vals1 = vals_onset2;
                features2 = (vals1 - stats1.stat_trials_mean_mean)./stats1.stat_trials_mean_sem;
            elseif strcmpi(feature, 'onset resp thresh')
                sel_cells2 = resp_cells_onset1;
                features2 = cat(1, stats1.onset_resp_thresh);
            elseif strcmpi(feature, 'offset mag')
                sel_cells2 = resp_cells_offset1;
                features2 = vals_offset2;
            elseif strcmpi(feature, 'offset mag z')
                sel_cells2 = resp_cells_offset1;
                vals1 = vals_offset2;
                features2 = (vals1 - stats1.stat_trials_mean_mean)./stats1.stat_trials_mean_sem;
            elseif strcmpi(feature, 'offset resp thresh')
                sel_cells2 = resp_cells_offset1;
                features2 = cat(1, stats1.offset_resp_thresh);
            end
            sel_cells{n_dset, n_gr, n_tn} = sel_cells2;
            features{n_dset, n_gr, n_tn} = features2;
            for n_reg = 1:num_reg
                reg_idx = logical(sum(area_labels0 == region_num(n_reg,:),2));
                cell_idx = logical(reg_idx .* sel_cells2);
                features_proc{n_dset, n_gr, n_reg, n_tn} = features2(cell_idx);
            end
        end
    end
end

end