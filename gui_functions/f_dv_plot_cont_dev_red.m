function f_dv_plot_cont_dev_red(app)

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupButtonGroup.SelectedObject.Text, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end

num_dsets = size(data,1);

tn_all2 = [18 29 20;28 19 30];

data_all = cell(num_dsets,1);
for n_dset = 1:num_dsets
    data1 = cell(2,1);
    for n_fl = 1:2
        stats1 = cat(1,data(n_dset,:).stats{n_pl});
        
        peak_val_all = stats1.peak_val_all(:, tn_all2(n_fl,:));
        
        resp_cells_all = f_dv_get_resp_cells(stats1, tn_all2(n_fl,:), app.LimitresptrialsCheckBox.Value, app.RespthreshEditField.Value);
        resp_cells = logical(sum(cat(2,resp_cells_all{:}),2));
        
        data1{n_fl} = peak_val_all(resp_cells,:);
    end
    data_all{n_dset} = cat(1,data1{:});
end

data_all2 = cat(1,data_all{:});


max_pt = max(data_all2(:));
min_pt = round(max_pt*.1,2);

f1 = figure;
plot3(data_all2(:,1), data_all2(:,2), data_all2(:,3), '.')
xlabel('cont'); ylabel('red'); zlabel('dev');
f1.Children.XDir = 'reverse';
f1.Children.YDir = 'reverse';
f1.Children.XLim = [-min_pt max_pt+min_pt];
f1.Children.YLim = [-min_pt max_pt+min_pt];
f1.Children.ZLim = [-min_pt max_pt+min_pt];
title(sprintf('%s peaks cont vs dev vs red', title_tag))

f2 = figure;
plot(data_all2(:,1), data_all2(:,3), '.');
xlabel('cont'); ylabel('dev');
f2.Children.XLim = [-min_pt max_pt+min_pt];
f2.Children.YLim = [-min_pt max_pt+min_pt];
title(sprintf('%s peaks cont vs dev', title_tag))

f3 = figure;
plot(data_all2(:,2), data_all2(:,3), '.');
xlabel('red'); ylabel('dev');
f3.Children.XLim = [-min_pt max_pt+min_pt];
f3.Children.YLim = [-min_pt max_pt+min_pt];
title(sprintf('%s peaks red vs dev', title_tag))

f4 = figure;
plot(data_all2(:,2), data_all2(:,1), '.');
xlabel('red'); ylabel('cont');
f4.Children.XLim = [-min_pt max_pt+min_pt];
f4.Children.YLim = [-min_pt max_pt+min_pt];
title(sprintf('%s peaks red vs cont', title_tag))

end