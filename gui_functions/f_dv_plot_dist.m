function f_dv_plot_dist(app)

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

num_dsets = numel(data.experiment);

tn_all = f_dv_get_trial_number(app);
[num_gr, num_tn] = size(tn_all);

features2 = cell(2,1);
sel_cells2 = cell(2,1);
for n_gr = 1:num_gr
    [features2{n_gr}, sel_cells2{n_gr}] = f_dv_get_feature(app, app.plotfeatureDropDown.Value, data, tn_all(n_gr,:), n_pl);
end

features1 = cat(1, features2{:});
sel_cells1 = cat(1, sel_cells2{:});
tn1 = tn_all(1,:);
% features1 = cell(1, num_tn);
% resp_cells1 = cell(1, num_tn);
% for n_tn = 1:num_tn
%     features1{n_tn} = cat(1,features3{:,n_tn});
%     resp_cells1{n_tn} = cat(1,resp_cells3{:,n_tn});
% end

plot_lims = f_str_to_array(app.plot_BaserespwinEditField.Value);
n_bins = round(diff(plot_lims)/data.cdata{1}.volume_period*1000/2);
figure; hold on; axis tight;
if app.MarginalizedistCheckBox.Value
    features_pool2 = cat(1, features1{:});
    resp_cells_pool2 = cat(1, sel_cells1{:});
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
    feat_pool = cell(1, num_tn);
    lab_pool = cell(1, num_tn);
    for n_tn = 1:num_tn
        features2 = cat(1,features1{:, n_tn});
        sel_cells2 = cat(1,sel_cells1{:, n_tn});
        feat_pool{n_tn} = features2(sel_cells2);
        lab_pool{n_tn} = repmat(tn1(n_tn), sum(sel_cells2),1);
        if sum(sel_cells1{n_tn})
            color2 = app.ops.context_types_all_colors2{tn1(n_tn)};
            if strcmpi(app.plottypeDropDown.Value, 'kde')
                [f, xi] = ksdensity(features2(sel_cells2), 'Function', 'pdf');
                plot(xi, f/sum(f), 'color', color2, 'LineWidth', 2);
            elseif strcmpi(app.plottypeDropDown.Value, 'ecdf')
                [f, xi] = ecdf(features2(sel_cells2));
                plot(xi, f, 'color', color2, 'LineWidth', 2);
            elseif strcmpi(app.plottypeDropDown.Value, 'histogram')
                histogram(features2(sel_cells2), linspace(plot_lims(1),plot_lims(2),n_bins), 'Normalization', 'probability');
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
    ylim([0 1.05])
end
title_tag2 = sprintf('%s, %s, %s, %s', title_tag, app.plotfeatureDropDown.Value, app.plottypeDropDown.Value, app.ResponsivecellsselectDropDown.Value);
title(title_tag2, 'interpreter', 'none');
legend(app.ops.context_types_labels_trim{tn1});

if ~app.MarginalizedistCheckBox.Value
    [p_all, tbl_all, stats_all]  = anova1(cat(1, feat_pool{:}),cat(1, lab_pool{:}), 'off');
    title_tag4 = sprintf('%s; stats', title_tag2);
    f_dv_plot_anova1(p_all, tbl_all, stats_all, title_tag4);
end
end