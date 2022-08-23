function tn = f_dv_get_trial_number(app)

if strcmpi(app.trialtypeDropDown.Value, 'all')
    tn = [1:10 19 20 29 30]; % 18 28
elseif strcmpi(app.trialtypeDropDown.Value, 'Freqs')
    tn = 1:10;
elseif strcmpi(app.trialtypeDropDown.Value, 'Context')
    tn = [18 29 20];
elseif strcmpi(app.trialtypeDropDown.Value, 'Context_flip')
    tn = [28 19 30];
elseif strcmpi(app.trialtypeDropDown.Value, 'Context_both')
    tn = [18 29 20 28 19 30];
elseif strcmpi(app.trialtypeDropDown.Value, 'Cont_dev')
    tn = [18 20];
elseif strcmpi(app.trialtypeDropDown.Value, 'Cont_dev_flip')
    tn = [28 30];
else
    tn = find(strcmpi(app.ops.context_types_labels, app.trialtypeDropDown.Value));
end 

end