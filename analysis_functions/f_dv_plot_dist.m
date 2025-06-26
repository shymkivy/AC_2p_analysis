function f_dv_plot_dist(app)

% app stuff
ops = app.ops;
params = f_dv_gather_params(app);

tn_all = f_dv_get_trial_number(params);
[data, title_tag] = f_dv_get_data_by_mouse_selection(app.data, params);
[region_num, reg_tag, leg_list] = f_dv_get_region_sel_val(params, ops);

[num_gr, num_tn] = size(tn_all);

[features2, sel_cells2] = f_dv_get_feature(params.plot_feature, tn_all, data, params, ops);
features1 = reshape(features2, [], num_tn);
sel_cells1 = reshape(sel_cells2, [], num_tn);

tn1 = tn_all(1,:);
% features1 = cell(1, num_tn);
% resp_cells1 = cell(1, num_tn);
% for n_tn = 1:num_tn
%     features1{n_tn} = cat(1,features3{:,n_tn});
%     resp_cells1{n_tn} = cat(1,resp_cells3{:,n_tn});
% end

num_bins = ceil(diff(params.plot_lims)/data.cdata{1}.volume_period*1000);
bin_locs = linspace(params.plot_lims(1),params.plot_lims(2),num_bins+1);
sm_fac = params.kde_smooth_factor;

figure; hold on; axis tight;
if sum(strcmpi(params.plot_type, {'kde', 'histogram', 'hist-kde'}))
    if params.plot_stim
        plot_times = 0:1:params.plot_lims(2);
        for n_pl = 1:numel(plot_times)
            r1 = rectangle('Position', [plot_times(n_pl) 0 0.5 1]);
            if strcmpi(ops.experiment_type, 'missmatch_grating')
                r1.FaceColor = [ops.context_types_all_colors2{params.stim_freq_color} params.stim_transparancy]; % +(n_pl-1)*2
                r1.EdgeColor = [ops.context_types_all_colors2{params.stim_freq_color} params.stim_transparancy]; % +(n_pl-1)*2
            else
                r1.FaceColor = [ops.context_types_all_colors2{params.stim_freq_color+(n_pl-1)*2} params.stim_transparancy]; 
                r1.EdgeColor = [ops.context_types_all_colors2{params.stim_freq_color+(n_pl-1)*2} params.stim_transparancy]; 
            end
        end
    end
end
ymax = 0;

if params.marginalize_dist
    features_pool2 = cat(1, features1{:});
    resp_cells_pool2 = cat(1, sel_cells1{:});
    data_feat = features_pool2(resp_cells_pool2);
    if sum(resp_cells_pool2)
        leg1 = cell(1,1);
        if strcmpi(params.plot_type, 'kde')
            [f1, xi] = if_kde_wrap(data_feat, bin_locs, sm_fac);
            leg1{1} = plot(xi, f1, 'LineWidth', 2, 'color', 'k');
            ymax = max([max(f1), ymax]);
            leg_lab = {'kde'};
        elseif strcmpi(params.plot_type, 'ecdf')
            [f, xi] = ecdf(data_feat);
            leg1{1} = plot(xi, f, 'LineWidth', 2, 'color', 'k');
            ymax = 1/1.1;
            leg_lab = {'ecdf'};
        elseif strcmpi(params.plot_type, 'histogram')
            h1 = histogram(data_feat, bin_locs, 'Normalization', 'probability');
            h1.FaceColor = [0.5 0.5 0.5];
            ymax = max([max(h1.Values), ymax]);
            leg1{1} = h1;
            leg_lab = {'hist'};
        elseif strcmpi(params.plot_type, 'hist-kde')
            leg1 = cell(2,1);
            h1 = histogram(data_feat, bin_locs, 'Normalization', 'probability'); % 
            h1.FaceColor = [0.5 0.5 0.5];
            leg1{1} = h1;
            [f1, xi] = if_kde_wrap(data_feat, bin_locs, sm_fac);
            leg1{2} = plot(xi, f1, 'LineWidth', 2, 'color', 'k');
            ymax = max([max(f1), max(h1.Values), ymax]);
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

        data_feat = features2(sel_cells2);

        if sum(sel_cells1{n_tn})
            has_data(n_tn) = 1;
            color2 = ops.context_types_all_colors2{tn1(n_tn)};
            if strcmpi(params.plot_type, 'kde')
                [f1, xi] = if_kde_wrap(data_feat, bin_locs, sm_fac);
                leg1{n_tn} = plot(xi, f1, 'color', color2, 'LineWidth', 2);
                ymax = max([max(f1), ymax]);
            elseif strcmpi(params.plot_type, 'ecdf')
                [f, xi] = ecdf(data_feat);
                leg1{n_tn} = plot(xi, f, 'color', color2, 'LineWidth', 2);
                ymax = 1/1.1;
            elseif strcmpi(params.plot_type, 'bar')
                mean1 = mean(data_feat);
                sem1 = std(data_feat)/sqrt(numel(data_feat)-1);
                leg1{n_tn} = bar(n_tn, mean1);
                leg1{n_tn}.FaceColor = color2;
                errorbar(n_tn, mean1, sem1, 'k');
                ymax = max([mean1+sem1, ymax]);
            elseif strcmpi(params.plot_type, 'histogram')
                h1 = histogram(data_feat, bin_locs, 'Normalization', 'probability');
                h1.FaceColor = color2;
                leg1{n_tn} = h1;
                ymax = max([max(h1.Values), ymax]);
            elseif strcmpi(params.plot_type, 'hist-kde')
                h1 = histogram(data_feat, bin_locs, 'Normalization', 'probability');
                h1.FaceColor = color2;
                leg1{n_tn} = h1;
                [f1, xi] = if_kde_wrap(data_feat, bin_locs, sm_fac);
                plot(xi, f1, 'color', color2, 'LineWidth', 2);
                ymax = max([max(f1), max(h1.Values), ymax]);
            end
        end
    end
end
ylim ([0 ymax*1.1])
ylabel('Fraction')
if strcmpi(params.plot_feature, 'peak loc')
    xlim(params.plot_lims);
    xlabel('Time, sec');
else
    xlabel(params.plot_feature);
end
title_tag2 = sprintf('%s, %s, %s, %s; smf=%.1f', title_tag, params.plot_feature, params.plot_type, params.responsive_cells_select, sm_fac);
title(title_tag2, 'interpreter', 'none');
if params.marginalize_dist
    legend([leg1{:}], leg_lab{:});
else
    legend([leg1{has_data}], ops.context_types_labels_trim2{tn1(has_data)});
end

if ~params.marginalize_dist
    [p_all, tbl_all, stats_all]  = anova1(cat(1, feat_pool{:}),cat(1, lab_pool{:}), 'off');
    title_tag4 = sprintf('%s; stats', title_tag2);
    f_dv_plot_anova1(p_all, tbl_all, stats_all, title_tag4, ops.context_types_labels_trim2(tn1(has_data)));
end
end


function [f1, xi] = if_kde_wrap(data_feat, bin_locs, sm_fac)

if ~exist('sm_fac', 'var')
    sm_fac = 1;
end

num_bins = numel(bin_locs)-1;
[f, xi] = ksdensity(data_feat, 'Function', 'pdf', 'NumPoints', 200, 'Bandwidth', 1/(num_bins-2)*sm_fac); % 
h_data = histcounts(data_feat, bin_locs)/numel(data_feat);
%idx1 = and(xi>plot_lims(1), xi<plot_lims(2));
idx1 = and(xi>max([min(data_feat),bin_locs(1)]), xi<min([max(data_feat),bin_locs(end)]));
f1 = f/(sum(f(idx1))/sum(idx1))/sum(logical(h_data))*sum(h_data);

end
