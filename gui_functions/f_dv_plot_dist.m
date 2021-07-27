function f_dv_plot_dist(app)

n_pl = app.mplSpinner.Value;
[data, title_tag] = f_dv_get_data_by_mouse_selection(app);
num_dsets = numel(data.experiment);

tn_all = f_dv_get_trial_number(app);
num_tn = numel(tn_all);

features1 = cell(num_dsets, num_tn);
for n_dset = 1:num_dsets
    
    stats1 = data(n_dset,:).stats{n_pl};
    
    if app.ConverttoZCheckBox.Value
        pop_mean_val = stats1.pop_mean_val;
        pop_z_factor = stats1.pop_z_factor;
    else
        pop_mean_val = zeros(stats1.num_cells,1);
        pop_z_factor = ones(stats1.num_cells,1);
    end
    
    for n_tn = 1:num_tn
        tn1 = tn_all(n_tn);
        
        cell_is_resp = stats1.cell_is_resp(:,tn1);
        
        if app.RespthreshEditField.Value
            resp_mag = stats1.peak_val_all(:,tn1);
            resp_mag = (resp_mag - pop_mean_val)./pop_z_factor;
            cell_is_resp = logical(cell_is_resp.*(resp_mag > app.RespthreshEditField.Value));
        end
        
        if strcmpi(app.plotfeatureDropDown.Value, 'peak loc')
            features1{n_dset, n_tn} = stats1.peak_t_all(cell_is_resp,tn1);
        elseif strcmpi(app.plotfeatureDropDown.Value, 'resp mag')
            features1{n_dset, n_tn} = stats1.peak_val_all(cell_is_resp,tn1);
            if app.ConverttoZCheckBox.Value
                features1{n_dset, n_tn} = (features1{n_dset, n_tn} - pop_mean_val(cell_is_resp))./pop_z_factor(cell_is_resp);
            end
        end
    end
end

features_pool = cell(1, num_tn);
for n_tn = 1:num_tn
    features_pool{n_tn} = cat(1, features1{:,n_tn});
end

figure; hold on; axis tight;
if app.MarginalizedistCheckBox.Value
    features_pool2 = cat(1, features_pool{:});
    if numel(features_pool2)
        if strcmpi(app.plottypeDropDown.Value, 'kde')
            [f, xi] = ksdensity(features_pool2);
            plot(xi, f, 'LineWidth', 2);
        elseif strcmpi(app.plottypeDropDown.Value, 'ecdf')
            [f, xi] = ecdf(features_pool2);
            plot(xi, f, 'LineWidth', 2);
        elseif strcmpi(app.plottypeDropDown.Value, 'histogram')
            histogram(features_pool2);
        end
    end
else
    for n_tn = 1:num_tn
        if numel(features_pool{n_tn})
            color2 = app.ops.context_types_all_colors2{tn_all(n_tn)};
            if strcmpi(app.plottypeDropDown.Value, 'kde')
                [f, xi] = ksdensity(features_pool{n_tn});
                plot(xi, f, 'color', color2, 'LineWidth', 2);
            elseif strcmpi(app.plottypeDropDown.Value, 'ecdf')
                [f, xi] = ecdf(features_pool{n_tn});
                plot(xi, f, 'color', color2, 'LineWidth', 2);
            elseif strcmpi(app.plottypeDropDown.Value, 'histogram')
                histogram(features_pool{n_tn});
            end
        end
    end
end
title(title_tag, 'interpreter', 'none')


    
end