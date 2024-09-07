function [features, sel_cells] = f_dv_get_feature(app, feature, data, tn_all, n_pl)

num_dsets = size(data,1);
num_tn = numel(tn_all);

features = cell(num_dsets, num_tn);
sel_cells = cell(num_dsets, num_tn);
for n_dset = 1:num_dsets
    stats1 = cat(1,data(n_dset,:).stats{n_pl});
    % first get responsive cells, then extract features
    
    sel_cells_peak = f_dv_get_resp_vals_cells(app, stats1, tn_all, 'peaks');
    sel_cells_onset = f_dv_get_resp_vals_cells(app, stats1, tn_all, 'onset');
    sel_cells_offset = f_dv_get_resp_vals_cells(app, stats1, tn_all, 'offset');
    
    for n_tn = 1:num_tn
        tn1 = tn_all(n_tn);
        
        sel_cells_peak1 = sel_cells_peak(:, n_tn);
        sel_cells_onset1 = sel_cells_onset(:, n_tn);
        sel_cells_offset1 = sel_cells_offset(:, n_tn);
        
        if strcmpi(feature, 'peak resp mag')
            sel_cells{n_dset, n_tn} = sel_cells_peak1;
            vals1 = cat(1, stats1.peak_vals);
            features{n_dset, n_tn} = vals1(:,tn1);
        elseif strcmpi(feature, 'peak resp mag z')
            sel_cells{n_dset, n_tn} = sel_cells_peak1;
            vals1 = cat(1, stats1.peak_vals);
            vals2 = vals1(:,tn1);
            features{n_dset, n_tn} = (vals2 - stats1.stat_trials_mean_mean)./stats1.stat_trials_mean_sem;
        elseif strcmpi(feature, 'peak loc')
            sel_cells{n_dset, n_tn} = sel_cells_peak1;
            vals1 = cat(1, stats1.peak_loc);
            features{n_dset, n_tn} = vals1(:,tn1);
        elseif strcmpi(feature, 'peak resp thresh')
            sel_cells{n_dset, n_tn} = sel_cells_peak1;
            features{n_dset, n_tn} = cat(1, stats1.peak_resp_thresh);
        elseif strcmpi(feature, 'stat trials mean sem')
            sel_cells{n_dset, n_tn} = sel_cells_peak1;
            features{n_dset, n_tn} = cat(1, stats1.stat_trials_mean_sem);
        elseif strcmpi(feature, 'onset mag')
            sel_cells{n_dset, n_tn} = sel_cells_onset1;
            vals1 = cat(1, stats1.onset_vals);
            features{n_dset, n_tn} = vals1(:,tn1);
        elseif strcmpi(feature, 'onset mag z')
            sel_cells{n_dset, n_tn} = sel_cells_onset1;
            vals1 = cat(1, stats1.onset_vals);
            vals2 = vals1(:,tn1);
            features{n_dset, n_tn} = (vals2 - stats1.stat_trials_mean_mean)./stats1.stat_trials_mean_sem;
        elseif strcmpi(feature, 'onset resp thresh')
            sel_cells{n_dset, n_tn} = sel_cells_onset1;
            features{n_dset, n_tn} = cat(1, stats1.onset_resp_thresh);
        elseif strcmpi(feature, 'offset mag')
            sel_cells{n_dset, n_tn} = sel_cells_offset1;
            vals1 = cat(1, stats1.offset_vals);
            features{n_dset, n_tn} = vals1(:,tn1);
        elseif strcmpi(feature, 'offset mag z')
            sel_cells{n_dset, n_tn} = sel_cells_offset1;
            vals1 = cat(1, stats1.offset_vals);
            vals2 = vals1(:,tn1);
            features{n_dset, n_tn} = (vals2 - stats1.stat_trials_mean_mean)./stats1.stat_trials_mean_sem;
        elseif strcmpi(feature, 'offset resp thresh')
            sel_cells{n_dset, n_tn} = sel_cells_offset1;
            features{n_dset, n_tn} = cat(1, stats1.offset_resp_thresh);
        end
    end
end

end