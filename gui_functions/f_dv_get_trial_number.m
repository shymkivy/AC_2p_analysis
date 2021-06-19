function tn = f_dv_get_trial_number(app)

if strcmpi(app.trialtypeDropDown.Value, 'all')
    tn = [1:10 19 20 29 30]; % 18 28
elseif strcmpi(app.trialtypeDropDown.Value, 'Freqs')
    tn = 1:10;
elseif strcmpi(app.trialtypeDropDown.Value, 'Context')
    tn = [18 19 20 28 29 30];
else
    tn = find(strcmpi(app.ops.context_types_labels, app.trialtypeDropDown.Value));
end 

end