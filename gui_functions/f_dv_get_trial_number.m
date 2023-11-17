function tn_out = f_dv_get_trial_number(app, trial_type_input, mmn_freq)

if ~exist('trial_type_input', 'var') || isempty(trial_type_input)
    trial_type_input = app.trialtypeDropDown.Value;
end

num_dsets = 1;
per_dset = 0;

if exist('mmn_freq', 'var')
    if iscell(mmn_freq)
        num_dsets = size(mmn_freq,1);
    end
end

if strcmpi(trial_type_input, 'all')
    tn = [1:10 18 19 20 28 29 30]; % 18 28
elseif strcmpi(trial_type_input, 'all_no_cont')
    tn = [1:10 19 20 29 30]; % 18 28
elseif strcmpi(trial_type_input, 'all_orig_sequnce')
    tn = [1:10 11:17 20 21:27 30]; % 18 28
elseif strcmpi(trial_type_input, 'Freqs')
    tn = 1:10;
elseif strcmpi(trial_type_input, 'Freqs -1')
    tn = 2:9;
elseif strcmpi(trial_type_input, 'Freqs -2')
    tn = 3:8;
elseif strcmpi(trial_type_input, 'Context')
    tn = [18 19 20];
elseif strcmpi(trial_type_input, 'Context_flip')
    tn = [28 29 30];
elseif strcmpi(trial_type_input, 'Context_both')
    tn = [18 19 20 28 29 30];
elseif strcmpi(trial_type_input, 'Context_both_comb')
    tn = [18 19 20; 28 29 30];
elseif strcmpi(trial_type_input, 'Cont_dev')
    tn = [18 20];
elseif strcmpi(trial_type_input, 'Cont_dev_flip')
    tn = [28 30];
elseif strcmpi(trial_type_input, 'Cont_dev_both')
    tn = [18 20 28 30];
elseif strcmpi(trial_type_input, 'Cont_dev_both_comb')
    tn = [18 20; 28 30];
elseif strcmpi(trial_type_input, 'Context +1')
    tn = [19 20 -1 0 1];
    mmn_idx = [0 0 2 2 2];
    per_dset = 1;
elseif strcmpi(trial_type_input, 'Context_flip +1')
    tn = [29 30 -1 0 1];
    mmn_idx = [0 0 1 1 1];
    per_dset = 1;
elseif strcmpi(trial_type_input, 'Context_both +1')
    tn = [19 20 -1 0 1 29 30 -1 0 1];
    mmn_idx = [0 0 2 2 2 0 0 1 1 1];
    per_dset = 1;
elseif strcmpi(trial_type_input, 'Context_both_comb +1')
    tn = [19 20 -1 0 1; 29 30 -1 0 1];
    mmn_idx = [0 0 2 2 2; 0 0 1 1 1];
    per_dset = 1;
else
    tn = find(strcmpi(app.ops.context_types_labels, trial_type_input));
end 

if per_dset
    tn_out = cell(num_dsets,1);
    for n_dset = 1:num_dsets
        tn1 = tn;
        mmn_freq1 = mmn_freq{n_dset};
        tn1(mmn_idx == 1) = tn(mmn_idx == 1) + mmn_freq1(1);
        tn1(mmn_idx == 2) = tn(mmn_idx == 2) + mmn_freq1(2);
        tn_out{n_dset} = tn1;
    end
else
    tn_out = tn;
end


end