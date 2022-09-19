function f_dv_plot_xy(app)

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

tn_all = f_dv_get_trial_number(app);

[featuresX, resp_cellsX] = f_dv_get_feature(app.xvarDropDown.Value, data, tn_all, n_pl, app.LimitresptrialsCheckBox.Value, app.RespthreshEditField.Value);
[featuresY, resp_cellsY] = f_dv_get_feature(app.yvarDropDown.Value, data, tn_all, n_pl, app.LimitresptrialsCheckBox.Value, app.RespthreshEditField.Value);

figure; hold on;
if app.MarginalizedistCheckBox.Value
    featuresX2 = cat(1,featuresX{:});
    resp_cellsX2 = cat(1,resp_cellsX{:});
    featuresY2 = cat(1,featuresY{:});
    resp_cellsY2 = cat(1,resp_cellsY{:});
    resp_all = logical(resp_cellsX2+resp_cellsY2);
    % OR
    %resp_all = logical(resp_cellsX2.*resp_cellsY2);
    plot(featuresX2(resp_all), featuresY2(resp_all), '.', 'color', [0 0 0]);
else
    for n_tn = 1:numel(tn_all)
        featuresX2 = cat(1,featuresX{:, n_tn});
        resp_cellsX2 = cat(1,resp_cellsX{:, n_tn});
        featuresY2 = cat(1,featuresY{:, n_tn});
        resp_cellsY2 = cat(1,resp_cellsY{:, n_tn});
        resp_all = logical(resp_cellsX2+resp_cellsY2);
        % OR
        %resp_all = logical(resp_cellsX2.*resp_cellsY2);
        plot(featuresX2(resp_all), featuresY2(resp_all), '.', 'color', app.ops.context_types_all_colors2{tn_all(n_tn)});
    end
end
xlabel(app.xvarDropDown.Value, 'interpreter', 'none')
ylabel(app.yvarDropDown.Value, 'interpreter', 'none')
title(title_tag, 'interpreter', 'none')

end