function tn = f_dv_get_trial_number(app)

if strcmpi(app.trialtypeDropDown.Value, 'all')
    tn = [1:10 18 19 20 28 29 30]; % 18 28
elseif strcmpi(app.trialtypeDropDown.Value, 'all_orig_sequnce')
    tn = [1:10 11:17 20 21:27 30]; % 18 28
elseif strcmpi(app.trialtypeDropDown.Value, 'Freqs')
    tn = 1:10;
elseif strcmpi(app.trialtypeDropDown.Value, 'Freqs -1')
    tn = 2:9;
elseif strcmpi(app.trialtypeDropDown.Value, 'Freqs -2')
    tn = 3:8;
elseif strcmpi(app.trialtypeDropDown.Value, 'Context')
    tn = [18 19 20];
elseif strcmpi(app.trialtypeDropDown.Value, 'Context_flip')
    tn = [28 29 30];
elseif strcmpi(app.trialtypeDropDown.Value, 'Context_both')
    tn = [18 19 20 28 29 30];
elseif strcmpi(app.trialtypeDropDown.Value, 'Context_both_comb')
    tn = [18 19 20; 28 29 30];
elseif strcmpi(app.trialtypeDropDown.Value, 'Cont_dev')
    tn = [18 20];
elseif strcmpi(app.trialtypeDropDown.Value, 'Cont_dev_flip')
    tn = [28 30];
elseif strcmpi(app.trialtypeDropDown.Value, 'Cont_dev_both')
    tn = [28 30 28 30];
elseif strcmpi(app.trialtypeDropDown.Value, 'Cont_dev_both_comb')
    tn = [28 30; 28 30];
else
    tn = find(strcmpi(app.ops.context_types_labels, app.trialtypeDropDown.Value));
end 

end