function f_dv_plot_dist_by_area(app)

% app stuff
ops = app.ops;
params = f_dv_gather_params(app);

tn_all = f_dv_get_trial_number(params);
[data, title_tag] = f_dv_get_data_by_mouse_selection(app.data, params);

[num_gr, num_tn] = size(tn_all);

[features0, sel_cells0, area_labels0] = f_dv_get_feature(params.plot_feature, tn_all, data, params, ops);
features1 = reshape(features0, [], num_tn);
sel_cells1 = reshape(sel_cells0, [], num_tn);
area_labels1 = reshape(area_labels0, [], 1);

tn1 = tn_all(1,:);

reg_all = app.ops.regions_to_analyze;
num_reg = numel(reg_all);

area_labels_pool = cat(1, area_labels1{:});

num_bins = ceil(diff(params.plot_lims)/data.cdata{1}.volume_period*1000);
bin_locs = linspace(params.plot_lims(1),params.plot_lims(2),num_bins+1);
sm_fac = params.kde_smooth_factor;

figure; hold on; axis tight;
if ~strcmpi(params.plot_type, 'ecdf')
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

feat_pool = cell(num_reg, num_tn);
lab_pool = cell(num_reg, num_tn);
if params.marginalize_dist
    use_reg = false(num_reg, 1);
    leg_pl = cell(num_reg, 1);
    leg_full = reg_all;
else
    use_reg = false(num_reg, num_tn);
    leg_full = cell(num_reg, num_tn);
    leg_pl = cell(num_reg, num_tn);
end
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

    if params.marginalize_dist
        features4 = features3(:);
        sel_cells4 = sel_cells3(:);
        data_feat = features4(sel_cells4);
        feat_pool{n_reg} = data_feat;
        lab_pool{n_reg} = repmat(tn1(n_tn), sum(sel_cells4),1);
        
        if sum(sel_cells4)
            color2 = ops.cond_colors{n_reg};
            use_reg(n_reg) = 1;
            if strcmpi(params.plot_type, 'kde')
                [f1, xi] = if_kde_wrap(data_feat, bin_locs, sm_fac);
                leg_pl{n_reg} = plot(xi, f1, 'color', color2, 'LineWidth', 2);
                ymax = max([max(f1), ymax]);
            elseif strcmpi(params.plot_type, 'ecdf')
                [f, xi] = ecdf(data_feat);
                leg_pl{n_reg} = plot(xi, f, 'color', color2, 'LineWidth', 2);
                ymax = 1/1.1;
            elseif strcmpi(params.plot_type, 'histogram')
                h1 = histogram(data_feat, bin_locs, 'Normalization', 'probability');
                h1.FaceColor = color2;
                leg_pl{n_reg} = h1;
                ymax = max([max(h1.Values), ymax]);
            elseif strcmpi(params.plot_type, 'hist-kde')
                h1 = histogram(data_feat, bin_locs, 'Normalization', 'probability'); % 
                h1.FaceColor = color2;
                leg_pl{n_reg} = h1;
                [f1, xi] = if_kde_wrap(data_feat, bin_locs, sm_fac);
                plot(xi, f1, 'color', color2, 'LineWidth', 2);
                ymax = max([max(f1), max(h1.Values), ymax]);
            end
        end
    else
        for n_tn = 1:num_tn
            features4 = features3(:, n_tn);
            sel_cells4 = sel_cells3(:, n_tn);
            data_feat = features4(sel_cells4);
            feat_pool{n_reg, n_tn} = data_feat;
            lab_pool{n_reg, n_tn} = repmat(n_reg*100+tn1(n_tn), sum(sel_cells4),1);
            leg_full{n_reg, n_tn} = sprintf('%s %s', reg_all{n_reg}, ops.context_types_labels_trim2{tn1(n_tn)});
            if sum(sel_cells4)
                use_reg(n_reg, n_tn) = 1;
                color2 = ops.context_types_all_colors2{tn1(n_tn)};
                linestyles2 = ops.cond_line_styles{n_reg};
                if strcmpi(params.plot_type, 'kde')
                    [f1, xi] = if_kde_wrap(data_feat, bin_locs, sm_fac);
                    leg_pl{n_reg, n_tn} = plot(xi, f1, 'color', color2, 'LineWidth', 2, 'LineStyle', linestyles2);
                    ymax = max([max(f1), ymax]);
                elseif strcmpi(params.plot_type, 'ecdf')
                    [f, xi] = ecdf(data_feat);
                    leg_pl{n_reg, n_tn} = plot(xi, f, 'color', color2, 'LineWidth', 2, 'LineStyle', linestyles2);
                    ymax = 1/1.1;
                elseif strcmpi(params.plot_type, 'histogram')
                    h1 = histogram(data_feat, bin_locs, 'Normalization', 'probability');
                    h1.FaceColor = color2;
                    leg_pl{n_reg, n_tn} = h1;
                elseif strcmpi(params.plot_type, 'hist-kde')
                    h1 = histogram(data_feat, bin_locs, 'Normalization', 'probability'); % 
                    h1.FaceColor = color2;
                    [f1, xi] = if_kde_wrap(data_feat, bin_locs, sm_fac);
                    leg_pl{n_reg, n_tn} = plot(xi, f1, 'color', color2, 'LineWidth', 2, 'LineStyle', linestyles2);
                    ymax = max([max(f1), max(h1.Values), ymax]);
                end
            end
        end
    end
end
legend_all = leg_full(use_reg);
leg_pl2 = leg_pl(use_reg);
ylim ([0 ymax*1.1])
ylabel('Fraction')
if strcmpi(params.plot_feature, 'peak loc')
    xlim(params.plot_lims);
    xlabel('Time, sec');
else
    xlabel(params.plot_feature);
end
title_tag2 = sprintf('%s, %s, %s; smf=%.1f', params.plot_feature, params.plot_type, params.responsive_cells_select, sm_fac);
title(title_tag2, 'interpreter', 'none');
legend([leg_pl2{:}], legend_all);

if ~params.marginalize_dist
    [p_all, tbl_all, stats_all]  = anova1(cat(1, feat_pool{:}),cat(1, lab_pool{:}), 'off');
    title_tag4 = sprintf('%s; stats', title_tag2);
    f_dv_plot_anova1(p_all, tbl_all, stats_all, title_tag4, legend_all);
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
