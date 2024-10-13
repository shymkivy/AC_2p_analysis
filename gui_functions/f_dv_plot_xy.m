function f_dv_plot_xy(app)

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

tn_all = f_dv_get_trial_number(app);

[featuresX, sel_cellsX] = f_dv_get_feature(app, app.xvarDropDown.Value, tn_all);
[featuresY, sel_cellsY] = f_dv_get_feature(app, app.yvarDropDown.Value, tn_all);

title_tagst = '';
f1 = figure; hold on; axis square;
if app.MarginalizedistCheckBox.Value
    featuresX2 = cat(1,featuresX{:});
    sel_cellsX2 = cat(1,sel_cellsX{:});
    featuresY2 = cat(1,featuresY{:});
    sel_cellsY2 = cat(1,sel_cellsY{:});
    sel_all = logical(sel_cellsX2+sel_cellsY2);
    % OR
    %resp_all = logical(resp_cellsX2.*resp_cellsY2);
    
    plot(featuresX2(sel_all), featuresY2(sel_all), '.', 'color', [0 0 0]);

    if app.plotstatsCheckBox.Value
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
        plot(featuresX2(sel_all), featuresY2(sel_all), '.', 'color', app.ops.context_types_all_colors2{tn_all(n_tn)});
        num_pts = num_pts + sum(sel_all);
    end
end
xlim1 = f1.Children.XLim;

% x_lab = xlim1(1):0.1:xlim1(2);
% plot(x_lab, x_lab, '--r');

xlabel(app.xvarDropDown.Value, 'interpreter', 'none')
ylabel(app.yvarDropDown.Value, 'interpreter', 'none')
title_tag2 = sprintf('%s; %s; %s; %dpts; %s', title_tag, app.trialtypeDropDown.Value, app.ResponsivecellsselectDropDown.Value, num_pts, title_tagst);
title(title_tag2, 'interpreter', 'none')


end