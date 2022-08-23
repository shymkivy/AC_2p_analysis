function resp_cells = f_dv_get_resp_cells(stats1, tn_all, resp_thr, limit_resp_trials)

resp_thr2 = stats1.resp_thresh(:,1)*resp_thr;

num_trials = numel(tn_all);

resp_cells = cell(num_trials,1);
for n_tn = 1:num_trials
    tn1 = tn_all(n_tn);
    if limit_resp_trials
        peak_vals = stats1.peak_val_all(:,tn1);
        resp_cells{n_tn} = logical(sum(peak_vals > resp_thr2,2));
    else
        peak_vals = stats1.peak_val_all(:,tn_all);
        resp_cells{n_tn} = logical(sum(peak_vals > resp_thr2,2));
    end
end

end