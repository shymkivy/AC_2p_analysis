function f_dv_plot_dist_by_area(app)

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

tn_all = f_dv_get_trial_number(app);
num_tn = numel(tn_all);

[features1, resp_cells, area_labels] = f_dv_get_feature2(app, app.plotfeatureDropDown.Value, data, tn_all, n_pl);

reg_all = app.ops.regions_to_analyze;
num_reg = numel(reg_all);

area_labels_pool = cat(1, area_labels{:});

legend_all = {};

plot_lims = f_str_to_array(app.plot_BaserespwinEditField.Value);
n_bins = round(diff(plot_lims)/data.cdata{1}.volume_period*1000/2);
figure; hold on; axis tight;
for n_reg = 1:num_reg
    area_idx = area_labels_pool == n_reg;
    num_area_cells = sum(area_idx);
    features3 = zeros(num_area_cells, num_tn);
    resp_cells3 = false(num_area_cells, num_tn);
    for n_tn = 1:num_tn
        features2 = cat(1,features1{:, n_tn});
        resp_cells2 = cat(1,resp_cells{:, n_tn});
        features3(:, n_tn) = features2(area_idx);
        resp_cells3(:, n_tn) = resp_cells2(area_idx);
    end

    if app.MarginalizedistCheckBox.Value
        features4 = features3(:);
        resp_cells4 = resp_cells3(:);
        if sum(resp_cells4)
            color2 = app.ops.cond_colors{n_reg};
            legend_all = [legend_all, {reg_all{n_reg}}];
            if strcmpi(app.plottypeDropDown.Value, 'kde')
                [f, xi] = ksdensity(features4(resp_cells4), 'Function', 'pdf');
                plot(xi, f/sum(f), 'color', color2, 'LineWidth', 2);
            elseif strcmpi(app.plottypeDropDown.Value, 'ecdf')
                [f, xi] = ecdf(features4(resp_cells4));
                plot(xi, f, 'color', color2, 'LineWidth', 2);
            elseif strcmpi(app.plottypeDropDown.Value, 'histogram')
                histogram(features4(resp_cells4), linspace(plot_lims(1),plot_lims(2),n_bins), 'Normalization', 'probability');
            end
        end
    else
        for n_tn = 1:num_tn
            features4 = features3(:, n_tn);
            resp_cells4 = resp_cells3(:, n_tn);
            if sum(resp_cells4)
                legend_all = [legend_all, {[reg_all{n_reg} ' ' app.ops.context_types_labels_trim{tn_all(n_tn)}]}];
                color2 = app.ops.context_types_all_colors2{tn_all(n_tn)};
                linestyles2 = app.ops.cond_line_styles{n_reg};
                if strcmpi(app.plottypeDropDown.Value, 'kde')
                    [f, xi] = ksdensity(features4(resp_cells4), 'Function', 'pdf');
                    plot(xi, f/sum(f), 'color', color2, 'LineWidth', 2, 'LineStyle', linestyles2);
                elseif strcmpi(app.plottypeDropDown.Value, 'ecdf')
                    [f, xi] = ecdf(features4(resp_cells4));
                    plot(xi, f, 'color', color2, 'LineWidth', 2, 'LineStyle', linestyles2);
                elseif strcmpi(app.plottypeDropDown.Value, 'histogram')
                    histogram(features2(resp_cells2), linspace(plot_lims(1),plot_lims(2),n_bins), 'Normalization', 'probability');
                end
            end
        end
    end
end

ylabel('Fraction')
if strcmpi(app.plotfeatureDropDown.Value, 'peak loc')
    xlim(plot_lims);
    xlabel('Time, sec');
else
    xlabel(app.plotfeatureDropDown.Value);
end
title(sprintf('%s, %s, %s, %s',title_tag, app.plotfeatureDropDown.Value, app.plottypeDropDown.Value, app.ResponsivecellsselectDropDown.Value), 'interpreter', 'none')
legend(legend_all)
end