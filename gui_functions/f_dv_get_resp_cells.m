function resp_cells = f_dv_get_resp_cells(stats1, tn_all, feature_type, limit_resp_trials, resp_thr)

if ~exist('feature_type', 'var')
    feature_type = 'peak';
end

resp_thr2 = stats1.resp_thresh_peak(:,1)*resp_thr;
num_trials = numel(tn_all);
resp_cells = cell(num_trials,1);

if strcmpi(feature_type, 'Peaks')
    vals1 = stats1.peak_val_all.*stats1.peak_in_resp_win;
elseif strcmpi(feature_type, 'Onset')
    vals1 = stats1.onset_vals;
elseif strcmpi(feature_type, 'Offset')
    vals1 = stats1.offset_vals;
%elseif strcmpi(feature_type, 'OnOff')
    
else
    error('feature type undefined');
end

for n_tn = 1:num_trials
    tn1 = tn_all(n_tn);
    if limit_resp_trials
        resp_cells{n_tn} = logical(sum(vals1(:,tn1) > resp_thr2,2));
    else
        resp_cells{n_tn} = logical(sum(vals1(:,tn_all) > resp_thr2,2));
    end
end

end