function [tn_out, con_idx] = f_dv_get_trial_number(app, mmn_freq)

trial_type_input = app.trialtypeDropDown.Value;

rel_freq = 0;

if strcmpi(trial_type_input, 'all')
    tn = [1:10 18 19 20 28 29 30]; % 18 28
    con_idx = {[1, 2, 3], [12, 11, 13], [14, 16, 15]};
elseif strcmpi(trial_type_input, 'all_no_cont')
    tn = [1:10 19 20 29 30]; % 18 28
    con_idx = {1:10};
elseif strcmpi(trial_type_input, 'all_orig_sequnce')
    tn = [1:10 11:17 20 21:27 30]; % 18 28
    con_idx = {1:10};
elseif strcmpi(trial_type_input, 'Freqs')
    tn = 1:10;
    con_idx = {1:10};
elseif strcmpi(trial_type_input, 'Freqs -1')
    tn = 2:9;
    con_idx = {1:8};
elseif strcmpi(trial_type_input, 'Freqs -2')
    tn = 3:8;
    con_idx = {1:6};
elseif strcmpi(trial_type_input, 'Context')
    tn = [18 19 20];
    con_idx = {[2, 1, 3]};
elseif strcmpi(trial_type_input, 'Context_flip')
    tn = [28 29 30];
    con_idx = {[2, 1, 3]};
elseif strcmpi(trial_type_input, 'Context_both')
    tn = [18 19 20 28 29 30];
    con_idx = {[2, 1, 3], [5, 4, 6]};
elseif strcmpi(trial_type_input, 'Context_both_comb')
    tn = [18 19 20; 28 29 30];
    con_idx = {[2, 1, 3]};
elseif strcmpi(trial_type_input, 'Cont_dev')
    tn = [18 20];
    con_idx = {[1, 2]};
elseif strcmpi(trial_type_input, 'Cont_dev_flip')
    tn = [28 30];
    con_idx = {[1, 2]};
elseif strcmpi(trial_type_input, 'Cont_dev_both')
    tn = [18 20 28 30];
    con_idx = {[1, 2], [3, 4]};
elseif strcmpi(trial_type_input, 'Cont_dev_both_comb')
    tn = [18 20; 28 30];
    con_idx = {[1, 2]};
elseif strcmpi(trial_type_input, 'Context +1')
    tn = [19 20 -1 0 1];
    mmn_idx = [0 0 2 2 2];
    rel_freq = 1;
    con_idx = {3:5, [1, 4, 2]};
elseif strcmpi(trial_type_input, 'Context_flip +1')
    tn = [29 30 -1 0 1];
    mmn_idx = [0 0 1 1 1];
    rel_freq = 1;
    con_idx = {3:5, [1, 4, 2]};
elseif strcmpi(trial_type_input, 'Context_both +1')
    tn = [19 20 -1 0 1 29 30 -1 0 1];
    mmn_idx = [0 0 2 2 2 0 0 1 1 1];
    rel_freq = 1;
    con_idx = {3:5, [1, 4, 2], 8:10, [6, 9, 7]};
elseif strcmpi(trial_type_input, 'Context_both_comb +1')
    tn = [19 20 -1 0 1; 29 30 -1 0 1];
    mmn_idx = [0 0 2 2 2; 0 0 1 1 1];
    rel_freq = 1;
    con_idx = {3:5, [1, 4, 2]};
elseif strcmpi(trial_type_input, 'Context +2')
    tn = [19 20 -2 -1 0 1 2];
    mmn_idx = [0 0 2 2 2 2 2];
    rel_freq = 1;
    con_idx = {3:7, [1, 5, 2]};
elseif strcmpi(trial_type_input, 'Context_flip +2')
    tn = [29 30 -2 -1 0 1 2];
    mmn_idx = [0 0 1 1 1 1 1];
    rel_freq = 1;
    con_idx = {3:7, [1, 5, 2]};
elseif strcmpi(trial_type_input, 'Context_both +2')
    tn = [19 20 -2 -1 0 1 2 29 30 -2 -1 0 1 2];
    mmn_idx = [0 0 2 2 2 2 2 0 0 1 1 1 1 1];
    rel_freq = 1;
    con_idx = {3:7, [1, 5, 2], 10:14, [8, 12, 9]};
elseif strcmpi(trial_type_input, 'Context_both_comb +2')
    tn = [19 20 -2 -1 0 1 2; 29 30 -2 -1 0 1 2];
    mmn_idx = [0 0 2 2 2 2 2; 0 0 1 1 1 1 1];
    rel_freq = 1;
    con_idx = {3:7, [1, 5, 2]};
else
    tn = find(strcmpi(app.ops.context_types_labels, trial_type_input));
end 

tn_out = tn;
if exist('mmn_freq', 'var')
    if rel_freq
        mmn1_idx = mmn_idx == 1;
        mmn2_idx = mmn_idx == 2;
    
        tn_out(mmn1_idx) = tn(mmn1_idx) + mmn_freq(1);
        tn_out(mmn2_idx) = tn(mmn2_idx) + mmn_freq(2);

        if sum(tn_out(mmn1_idx)>10)
            idx1 = logical((tn_out > 10) .* mmn1_idx);
            tn_out(idx1) = 0;
        elseif sum(tn_out(mmn2_idx)>10)
            idx1 = logical((tn_out > 10) .* mmn2_idx);
            tn_out(idx1) = 0;
        elseif sum(tn_out(mmn1_idx)<1)
            idx1 = logical((tn_out < 1) .* mmn1_idx);
            tn_out(idx1) = 0;
        elseif sum(tn_out(mmn2_idx)<1)
            idx1 = logical((tn_out < 1) .* mmn2_idx);
            tn_out(idx1) = 0;
        end
    end
end

end