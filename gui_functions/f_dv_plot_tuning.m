function f_dv_plot_tuning(app)

n_pl = app.mplSpinner.Value;
[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

tn_all = f_dv_get_trial_number(app);

num_dsets = numel(data.experiment);

resp_cell_all = zeros(num_dsets, numel(tn_all));
for n_tn = 1:numel(tn_all)
    for n_dset = 1:num_dsets
        tn1 = tn_all(n_tn);
        
        data1 = data(n_dset,:);

        resp_cells1 = data1.stats{1}{n_pl}.cell_is_resp(:,tn1);
        
        resp_cell_all(n_dset, n_tn) = sum(resp_cells1);
        
    end
end

categories = app.ops.context_types_labels(tn_all);

if app.poolgroupsCheckBox.Value
    resp1= sum(resp_cell_all,1);
else
    resp1 = resp_cell_all;
end

figure;
bar(categorical(categories,categories), resp1, 'EdgeColor',[219, 239, 255]/256,'LineWidth',1.5);
title(title_tag)

end