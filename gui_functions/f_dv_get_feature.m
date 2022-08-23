function feature_out = f_dv_get_feature(feature, data, tn_all, n_pl, limit_resp_trials, resp_thr)

num_dsets = size(data,1);
num_tn = numel(tn_all);

feature1 = cell(num_dsets, num_tn);
for n_dset = 1:num_dsets
    stats1 = data(n_dset,:).stats{n_pl};
    % first get responsive cells, then extract features
   
    resp_cells_all = f_dv_get_resp_cells(stats1, tn_all, limit_resp_trials, resp_thr);
    
    for n_tn = 1:num_tn
        tn1 = tn_all(n_tn);
        
        resp_cells = resp_cells_all{n_tn};
        
        if strcmpi(feature, 'resp mag')
            feature1{n_dset, n_tn} = stats1.peak_val_all(resp_cells,tn1);
        elseif strcmpi(feature, 'resp mag z')
            peak_vals = stats1.peak_val_all(resp_cells,tn1);
            feature1{n_dset, n_tn} = (peak_vals - stats1.trial_ave_val(resp_cells))./stats1.trial_sem_val(resp_cells);
        elseif strcmpi(feature, 'peak loc')
            feature1{n_dset, n_tn} = stats1.peak_t_all(resp_cells,tn1);
        elseif strcmpi(feature, 'resp thresh')
            feature1{n_dset, n_tn} = stats1.resp_thresh(resp_cells);
        elseif strcmpi(feature, 'trial sem val')
            feature1{n_dset, n_tn} = stats1.trial_sem_val(resp_cells);
        end
    end
end

feature_out = feature1;

end