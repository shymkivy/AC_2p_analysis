function f_dv_plot_xy(app)

% app stuff
ops = app.ops;
params = f_dv_gather_params(app);

tn_all = f_dv_get_trial_number(params);
[data, title_tag] = f_dv_get_data_by_mouse_selection(app.data, params);
%[region_num, reg_tag, leg_list] = f_dv_get_region_sel_val(params, ops);

[featuresX, sel_cellsX] = f_dv_get_feature(params.x_var, tn_all, data, params, ops);
[featuresY, sel_cellsY] = f_dv_get_feature(params.y_var, tn_all, data, params, ops);

title_tagst = '';
f1 = figure; hold on; axis square;
if params.marginalize_dist
    featuresX2 = cat(1,featuresX{:});
    sel_cellsX2 = cat(1,sel_cellsX{:});
    featuresY2 = cat(1,featuresY{:});
    sel_cellsY2 = cat(1,sel_cellsY{:});
    sel_all = logical(sel_cellsX2+sel_cellsY2);
    % OR
    %resp_all = logical(resp_cellsX2.*resp_cellsY2);
    
    plot(featuresX2(sel_all), featuresY2(sel_all), '.', 'color', [0 0 0]);

    if params.plot_stats
        [rho,pval] = corr(featuresX2(sel_all), featuresY2(sel_all));
        title_tagst = sprintf('rho=%.2f; p=%.2e', rho, pval);
    end

    num_pts = sum(sel_all);
else
    num_pts = 0;
    for n_tn = 1:numel(tn_all)
        featuresX2 = cat(1,featuresX{:, n_tn});
        sel_cellsX2 = cat(1,sel_cellsX{:, n_tn});
        featuresY2 = cat(1,featuresY{:, n_tn});
        sel_cellsY2 = cat(1,sel_cellsY{:, n_tn});
        sel_all = logical(sel_cellsX2+sel_cellsY2);
        % OR
        %resp_all = logical(resp_cellsX2.*resp_cellsY2);
        plot(featuresX2(sel_all), featuresY2(sel_all), '.', 'color', ops.context_types_all_colors2{tn_all(n_tn)});
        num_pts = num_pts + sum(sel_all);
    end
end
xlim1 = f1.Children.XLim;

% x_lab = xlim1(1):0.1:xlim1(2);
% plot(x_lab, x_lab, '--r');

xlabel(params.x_var, 'interpreter', 'none')
ylabel(params.y_var, 'interpreter', 'none')
title_tag2 = sprintf('%s; %s; %s; %dpts; %s', title_tag, params.trial_type, params.responsive_cells_select, num_pts, title_tagst);
title(title_tag2, 'interpreter', 'none')


end