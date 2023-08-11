function f_dv_plot_dist_by_area(app)

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

tn_all = f_dv_get_trial_number(app);
[num_gr, num_tn] = size(tn_all);

features0 = cell(num_gr,1);
sel_cells0 = cell(num_gr,1);
area_labels0 = cell(num_gr,1);


transp = app.StimtranspEditField.Value;
freq_col = app.stimcolorSpinner.Value;


for n_gr = 1:num_gr
    [features0{n_gr}, sel_cells0{n_gr}, area_labels0{n_gr}] = f_dv_get_feature2(app, app.plotfeatureDropDown.Value, data, tn_all(n_gr,:), n_pl);
end

features1 = cat(1, features0{:});
sel_cells1 = cat(1, sel_cells0{:});
area_labels1 = cat(1, area_labels0{:});
tn1 = tn_all(1,:);

reg_all = app.ops.regions_to_analyze;
num_reg = numel(reg_all);

area_labels_pool = cat(1, area_labels1{:});

legend_all = {};

plot_lims = f_str_to_array(app.plot_BaserespwinEditField.Value);
n_bins = ceil(diff(plot_lims)/data.cdata{1}.volume_period*1000);
figure; hold on; axis tight;
if app.PlotstimCheckBox.Value
    plot_times = 0:1:plot_lims(2);
    for n_pl = 1:numel(plot_times)
        r1 = rectangle('Position', [plot_times(n_pl) 0 0.5 1]);
        if strcmpi(app.ops.experiment_type, 'missmatch_grating')
            r1.FaceColor = [app.ops.context_types_all_colors2{freq_col} transp]; % +(n_pl-1)*2
            r1.EdgeColor = [app.ops.context_types_all_colors2{freq_col} transp]; % +(n_pl-1)*2
        else
            r1.FaceColor = [app.ops.context_types_all_colors2{freq_col+(n_pl-1)*2} transp]; 
            r1.EdgeColor = [app.ops.context_types_all_colors2{freq_col+(n_pl-1)*2} transp]; 
        end
    end
end
ymax = 0;

feat_pool = cell(num_reg, num_tn);
lab_pool = cell(num_reg, num_tn);
for n_reg = 1:num_reg
    area_idx = area_labels_pool == n_reg;
    num_area_cells = sum(area_idx);
    features3 = zeros(num_area_cells, num_tn);
    sel_cells3 = false(num_area_cells, num_tn);
    for n_tn = 1:num_tn
        features2 = cat(1,features1{:, n_tn});
        sel_cells2 = cat(1,sel_cells1{:, n_tn});
        features3(:, n_tn) = features2(area_idx);
        sel_cells3(:, n_tn) = sel_cells2(area_idx);
    end

    if app.MarginalizedistCheckBox.Value
        features4 = features3(:);
        sel_cells4 = sel_cells3(:);
        feat_pool{n_reg} = features4(sel_cells4);
        lab_pool{n_reg} = repmat(tn1(n_tn), sum(sel_cells4),1);
        if sum(sel_cells4)
            color2 = app.ops.cond_colors{n_reg};
            legend_all = [legend_all, {reg_all{n_reg}}];
            if strcmpi(app.plottypeDropDown.Value, 'kde')
                [f, xi] = ksdensity(features4(sel_cells4), 'Function', 'pdf', 'NumPoints', 200);
                idx1 = and(xi>plot_lims(1), xi<plot_lims(2));
                y1 = f/sum(f(idx1))/n_bins*sum(idx1);
                ymax = max([max(y1), ymax]);
                plot(xi, y1, 'color', color2, 'LineWidth', 2);
            elseif strcmpi(app.plottypeDropDown.Value, 'ecdf')
                [f, xi] = ecdf(features4(sel_cells4));
                plot(xi, f, 'color', color2, 'LineWidth', 2);
                ymax = 1;
            elseif strcmpi(app.plottypeDropDown.Value, 'histogram')
                h1 = histogram(features4(sel_cells4), linspace(plot_lims(1),plot_lims(2),n_bins), 'Normalization', 'probability');
                h1.FaceColor = color2;
                ymax = max([max(h1.Values), ymax]);
            end
        end
    else
        for n_tn = 1:num_tn
            features4 = features3(:, n_tn);
            sel_cells4 = sel_cells3(:, n_tn);
            feat_pool{n_reg, n_tn} = features4(sel_cells4);
            lab_pool{n_reg, n_tn} = repmat(n_reg*100+tn1(n_tn), sum(sel_cells4),1);
            if sum(sel_cells4)
                legend_all = [legend_all, {[reg_all{n_reg} ' ' app.ops.context_types_labels_trim{tn1(n_tn)}]}];
                color2 = app.ops.context_types_all_colors2{tn1(n_tn)};
                linestyles2 = app.ops.cond_line_styles{n_reg};
                if strcmpi(app.plottypeDropDown.Value, 'kde')
                    [f, xi] = ksdensity(features4(sel_cells4), 'Function', 'pdf', 'NumPoints', 200);
                    idx1 = and(xi>plot_lims(1), xi<plot_lims(2));
                    y1 = f/sum(f(idx1))/n_bins*sum(idx1);
                    plot(xi, y1, 'color', color2, 'LineWidth', 2, 'LineStyle', linestyles2);
                elseif strcmpi(app.plottypeDropDown.Value, 'ecdf')
                    [f, xi] = ecdf(features4(sel_cells4));
                    plot(xi, f, 'color', color2, 'LineWidth', 2, 'LineStyle', linestyles2);
                elseif strcmpi(app.plottypeDropDown.Value, 'histogram')
                    histogram(features2(sel_cells2), linspace(plot_lims(1),plot_lims(2),n_bins), 'Normalization', 'probability');
                end
            end
        end
    end
end
ylim ([0 ymax*1.1])
ylabel('Fraction')
if strcmpi(app.plotfeatureDropDown.Value, 'peak loc')
    xlim(plot_lims);
    xlabel('Time, sec');
else
    xlabel(app.plotfeatureDropDown.Value);
    ylim([0 1.05])
end
title_tag2 = sprintf('%s, %s, %s, %s', app.plotfeatureDropDown.Value, app.plottypeDropDown.Value, app.ResponsivecellsselectDropDown.Value);
title(title_tag2, 'interpreter', 'none');
legend(legend_all);

if ~app.MarginalizedistCheckBox.Value
    [p_all, tbl_all, stats_all]  = anova1(cat(1, feat_pool{:}),cat(1, lab_pool{:}), 'off');
    title_tag4 = sprintf('%s; stats', title_tag2);
    f_dv_plot_anova1(p_all, tbl_all, stats_all, title_tag4);
end

end