function f_dv_plot_dist(app)

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupButtonGroup.SelectedObject.Text, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

num_dsets = numel(data.experiment);

tn_all = f_dv_get_trial_number(app);
num_tn = numel(tn_all);

features1 = f_dv_get_feature(app.plotfeatureDropDown.Value, data, tn_all, n_pl, app.LimitresptrialsCheckBox.Value, app.RespthreshEditField.Value);

figure; hold on; axis tight;
if app.MarginalizedistCheckBox.Value
    features_pool2 = cat(1, features1{:});
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
        if numel(features1{n_tn})
            color2 = app.ops.context_types_all_colors2{tn_all(n_tn)};
            if strcmpi(app.plottypeDropDown.Value, 'kde')
                [f, xi] = ksdensity(features1{n_tn});
                plot(xi, f, 'color', color2, 'LineWidth', 2);
            elseif strcmpi(app.plottypeDropDown.Value, 'ecdf')
                [f, xi] = ecdf(features1{n_tn});
                plot(xi, f, 'color', color2, 'LineWidth', 2);
            elseif strcmpi(app.plottypeDropDown.Value, 'histogram')
                histogram(features1{n_tn});
            end
        end
    end
end
title(title_tag, 'interpreter', 'none')


    
end