function trial_types_ctx = f_dv_mark_tt_ctx(trial_types, MMN_freq, ops)

trial_types_ctx = zeros(size(trial_types,1),2);

if ~isempty(MMN_freq)
    % mark controls for ctx MMM2 = 170
    % idx1 = trial_types == MMN_freq(1);
    % trial_types_ctx(idx1,2) = 250;
    for n_freq = 1:10
        idx1 = trial_types == n_freq;
        trial_types_ctx(idx1,2) = n_freq - MMN_freq(1) + 250;
    end

    % idx1 = trial_types == MMN_freq(2);
    % trial_types_ctx(idx1,1) = 150;
    for n_freq = 1:10
        idx1 = trial_types == n_freq;
        trial_types_ctx(idx1,1) = n_freq - MMN_freq(2) + 150;
    end
end

redf_idx = logical(sum(trial_types == (200 + ops.redundent_pool_trials),2));
trial_types_ctx(redf_idx,1) = 260;

redf_idx = logical(sum(trial_types == (100 + ops.redundent_pool_trials),2));
trial_types_ctx(redf_idx,2) = 160;

end