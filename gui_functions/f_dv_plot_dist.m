function f_dv_plot_dist(app)

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

num_dsets = numel(data.experiment);

tn_all = f_dv_get_trial_number(app);
num_tn = numel(tn_all);

[features1, resp_cells] = f_dv_get_feature(app, app.plotfeatureDropDown.Value, data, tn_all, n_pl);

plot_lims = f_str_to_array(app.plot_BaserespwinEditField.Value);
n_bins = round(diff(plot_lims)/data.cdata{1}.volume_period*1000/2);
figure; hold on; axis tight;
if app.MarginalizedistCheckBox.Value
    features_pool2 = cat(1, features1{:});
    resp_cells_pool2 = cat(1, resp_cells{:});
    if sum(resp_cells_pool2)
        if strcmpi(app.plottypeDropDown.Value, 'kde')
            [f, xi] = ksdensity(features_pool2(resp_cells_pool2), 'Function', 'pdf');
            plot(xi, f/sum(f), 'LineWidth', 2);
        elseif strcmpi(app.plottypeDropDown.Value, 'ecdf')
            [f, xi] = ecdf(features_pool2(resp_cells_pool2));
            plot(xi, f, 'LineWidth', 2);
        elseif strcmpi(app.plottypeDropDown.Value, 'histogram')
            histogram(features_pool2(resp_cells_pool2), linspace(plot_lims(1),plot_lims(2),n_bins), 'Normalization', 'probability');
        end
    end
else
    for n_tn = 1:num_tn
        features2 = cat(1,features1{:, n_tn});
        resp_cells2 = cat(1,resp_cells{:, n_tn});
        if sum(resp_cells{n_tn})
            color2 = app.ops.context_types_all_colors2{tn_all(n_tn)};
            if strcmpi(app.plottypeDropDown.Value, 'kde')
                [f, xi] = ksdensity(features2(resp_cells2), 'Function', 'pdf');
                plot(xi, f/sum(f), 'color', color2, 'LineWidth', 2);
            elseif strcmpi(app.plottypeDropDown.Value, 'ecdf')
                [f, xi] = ecdf(features2(resp_cells2));
                plot(xi, f, 'color', color2, 'LineWidth', 2);
            elseif strcmpi(app.plottypeDropDown.Value, 'histogram')
                histogram(features2(resp_cells2), linspace(plot_lims(1),plot_lims(2),n_bins), 'Normalization', 'probability');
            end
        end
    end
end
ylabel('Fraction')
if strcmpi(app.plotfeatureDropDown.Value, 'peak loc')
    xlim(plot_lims);
    xlabel('Time, sec');
else
    xlabel('Magnitude');
end
title(sprintf('%s, %s, %s, %s',title_tag, app.plotfeatureDropDown.Value, app.plottypeDropDown.Value, app.ResponsivecellsselectDropDown.Value), 'interpreter', 'none')
    
end