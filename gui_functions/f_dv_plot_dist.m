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

transp = app.StimtranspEditField.Value;
freq_col = app.stimcolorSpinner.Value;

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
n_bins = ceil(diff(plot_lims)/data.cdata{1}.volume_period*1000);
figure; hold on; axis tight;

if ~strcmpi(app.plottypeDropDown.Value, 'ecdf')
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
end
ymax = 0;
if app.MarginalizedistCheckBox.Value
    features_pool2 = cat(1, features1{:});
    resp_cells_pool2 = cat(1, sel_cells1{:});
    if sum(resp_cells_pool2)
        if strcmpi(app.plottypeDropDown.Value, 'kde')
            leg1 = cell(1,1);
            [f, xi] = ksdensity(features_pool2(resp_cells_pool2), 'Function', 'pdf', 'NumPoints', 200);
            idx1 = and(xi>plot_lims(1), xi<plot_lims(2));
            y1 = f/sum(f(idx1))/n_bins*sum(idx1);
            leg1{1} = plot(xi, y1, 'LineWidth', 2, 'color', 'k');
            ymax = max([max(y1), ymax]);
            leg_lab = {'kde'};
        elseif strcmpi(app.plottypeDropDown.Value, 'ecdf')
            leg1 = cell(1,1);
            [f, xi] = ecdf(features_pool2(resp_cells_pool2));
            leg1{1} = plot(xi, f, 'LineWidth', 2, 'color', 'k');
            ymax = 1;
            leg_lab = {'ecdf'};
        elseif strcmpi(app.plottypeDropDown.Value, 'histogram')
            leg1 = cell(1,1);
            h1 = histogram(features_pool2(resp_cells_pool2), linspace(plot_lims(1),plot_lims(2),n_bins), 'Normalization', 'probability');
            h1.FaceColor = [0.5 0.5 0.5];
            ymax = max([max(h1.Values), ymax]);
            leg1{1} = h1;
            leg_lab = {'hist'};
        elseif strcmpi(app.plottypeDropDown.Value, 'hist-kde')
            leg1 = cell(2,1);
            h1 = histogram(features_pool2(resp_cells_pool2), linspace(plot_lims(1),plot_lims(2),n_bins), 'Normalization', 'probability');
            h1.FaceColor = [0.5 0.5 0.5];
            leg1{1} = h1;
            ymax = max([max(h1.Values), ymax]);
            [f, xi] = ksdensity(features_pool2(resp_cells_pool2), 'Function', 'pdf', 'NumPoints', 200, 'Bandwidth', 0.02);
            idx1 = and(xi>plot_lims(1), xi<plot_lims(2));
            y1 = f/sum(f(idx1))/n_bins*sum(idx1);
            leg1{2} = plot(xi, y1, 'LineWidth', 2, 'color', 'k');
            ymax = max([max(y1), ymax]);
            leg_lab = {'histogram', 'KDE'};
        end
    end
else
    has_data = false(num_tn,1);
    leg1 = cell(num_tn,1);
    feat_pool = cell(1, num_tn);
    lab_pool = cell(1, num_tn);
    for n_tn = 1:num_tn
        features2 = cat(1,features1{:, n_tn});
        sel_cells2 = cat(1,sel_cells1{:, n_tn});
        feat_pool{n_tn} = features2(sel_cells2);
        lab_pool{n_tn} = repmat(tn1(n_tn), sum(sel_cells2),1);
        if sum(sel_cells1{n_tn})
            has_data(n_tn) = 1;
            color2 = app.ops.context_types_all_colors2{tn1(n_tn)};
            if strcmpi(app.plottypeDropDown.Value, 'kde')
                [f, xi] = ksdensity(features2(sel_cells2), 'Function', 'pdf');
                idx1 = and(xi>plot_lims(1), xi<plot_lims(2));
                y1 = f/sum(f(idx1))/n_bins*sum(idx1);
                leg1{n_tn} = plot(xi, y1, 'color', color2, 'LineWidth', 2);
                ymax = max([max(y1), ymax]);
            elseif strcmpi(app.plottypeDropDown.Value, 'ecdf')
                [f, xi] = ecdf(features2(sel_cells2));
                leg1{n_tn} = plot(xi, f, 'color', color2, 'LineWidth', 2);
                ymax = 1;
            elseif strcmpi(app.plottypeDropDown.Value, 'histogram')
                h1 = histogram(features2(sel_cells2), linspace(plot_lims(1),plot_lims(2),n_bins), 'Normalization', 'probability');
                h1.FaceColor = color2;
                leg1{n_tn} = h1;
                ymax = max([max(h1.Values), ymax]);
            elseif strcmpi(app.plottypeDropDown.Value, 'hist-kde')
                h1 = histogram(features2(sel_cells2), linspace(plot_lims(1),plot_lims(2),n_bins), 'Normalization', 'probability');
                h1.FaceColor = color2;
                leg1{n_tn} = h1;
                [f, xi] = ksdensity(features2(sel_cells2), 'Function', 'pdf');
                idx1 = and(xi>plot_lims(1), xi<plot_lims(2));
                y1 = f/sum(f(idx1))/n_bins*sum(idx1);
                plot(xi, y1, 'color', color2, 'LineWidth', 2);
                ymax = max([max(y1), ymax]);
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
title_tag2 = sprintf('%s, %s, %s, %s', title_tag, app.plotfeatureDropDown.Value, app.plottypeDropDown.Value, app.ResponsivecellsselectDropDown.Value);
title(title_tag2, 'interpreter', 'none');
if app.MarginalizedistCheckBox.Value
    legend([leg1{:}], leg_lab{:});
else
    legend([leg1{has_data}], app.ops.context_types_labels_trim{tn1(has_data)});
end

if ~app.MarginalizedistCheckBox.Value
    [p_all, tbl_all, stats_all]  = anova1(cat(1, feat_pool{:}),cat(1, lab_pool{:}), 'off');
    title_tag4 = sprintf('%s; stats', title_tag2);
    f_dv_plot_anova1(p_all, tbl_all, stats_all, title_tag4);
end
end