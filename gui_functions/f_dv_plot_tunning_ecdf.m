function f_dv_plot_tunning_ecdf(app)
% can either for every cell get best tuning and use that, or for each cond
% get tuning


n_pl = app.mplSpinner.Value;

resp_cells = app.ddata.stats{1}{n_pl}.cell_is_resp;
peak_vals = app.ddata.stats{1}{n_pl}.peak_val_all;
%peak_vals(~resp_cells) = 0;

if app.ConverttoZCheckBox.Value
    pop_mean_val = app.ddata.stats{1}{n_pl}.pop_mean_val;
    pop_z_factor = app.ddata.stats{1}{n_pl}.pop_z_factor;
end

tn_all = f_dv_get_trial_number(app);

figure; hold on;
if app.CombinegroupsCheckBox.Value
    
    peak_vals2 = max(peak_vals(:,tn_all),[],2);
    if app.ConverttoZCheckBox.Value
        peak_vals2 = (peak_vals2 - pop_mean_val)./pop_z_factor;
    end
    resp_cells2 = logical(sum(resp_cells(:,tn_all),2));
    peak_vals2(~resp_cells2) = [];  
    peak_vals_sort = sort(peak_vals2);
    f = (1:numel(peak_vals_sort))/numel(peak_vals_sort);
    %[f, x] = ecdf(reliab_list);
    plot(peak_vals_sort,f, 'LineWidth', 2);
else
    for n_tn = 1:numel(tn_all)
        tn_temp = tn_all(n_tn); 
        peak_vals2 = peak_vals(:,tn_temp);
        if app.ConverttoZCheckBox.Value
            peak_vals2 = (peak_vals2 - pop_mean_val)./pop_z_factor;
        end
        
        resp_cells2 = resp_cells(:,tn_temp);
        peak_vals2(~resp_cells2) = [];        
        peak_vals_sort = sort(peak_vals2);
        f = (1:numel(peak_vals_sort))/numel(peak_vals_sort); 
        
        if tn_temp>10
            plot(peak_vals_sort,f, 'LineWidth', 2, 'Color', app.ops.context_types_all_colors2{tn_temp}, 'LineStyle', '--');
        else
            plot(peak_vals_sort,f, 'LineWidth', 2, 'Color', app.ops.context_types_all_colors2{tn_temp});
        end
    end
end

title(sprintf(' %s; Tuning mag ecdf; %s', app.ddata.experiment{1}), app.trialtypeDropDown.Value, 'Interpreter', 'none')
axis tight;
if app.ConverttoZCheckBox.Value
    xlabel('resp Z-score');
else
    xlabel('resp mag');
end

end