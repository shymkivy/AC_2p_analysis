function trial_numbers = f_tt_to_tn(trial_types, ops, convert_reds)

if ~exist('convert_reds', 'var') || isempty(convert_reds)
    convert_reds = 0;
end
if convert_reds
    for n_tr = 1:numel(trial_types)
        if (trial_types(n_tr) > 108) && (trial_types(n_tr) < 160)
            trial_types(n_tr) = 108;
        end
        if (trial_types(n_tr) > 208) && (trial_types(n_tr) < 260)
            trial_types(n_tr) = 208;
        end
    end
end

trial_numbers = zeros(numel(trial_types),1);
for n_tr = 1:numel(trial_types)
    tn = find(ops.context_types_all == trial_types(n_tr));
    if ~isempty(tn)
    	trial_numbers(n_tr) = tn;
    end
end
end