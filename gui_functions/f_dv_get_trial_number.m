function tn = f_dv_get_trial_number(app, trial_type_input)

if ~exist('man_input', 'var') || isempty(man_input)
    trial_type_input = app.trialtypeDropDown.Value;
end

if strcmpi(trial_type_input, 'all')
    tn = [1:10 18 19 20 28 29 30]; % 18 28
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
    tn = [28 30 28 30];
elseif strcmpi(trial_type_input, 'Cont_dev_both_comb')
    tn = [28 30; 28 30];
else
    tn = find(strcmpi(app.ops.context_types_labels, trial_type_input));
end 

end