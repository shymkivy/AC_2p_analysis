function f_dv_more_tuining_plot(app)

[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

if strcmpi(app.SelectdatagroupDropDown.Value, 'plane')
    n_pl = app.mplSpinner.Value;
else
    n_pl = 1:max([data.num_planes]);
end


stats1 = data.stats{1};


resp1 = stats1.peak_vals(:,1:10);



figure; imagesc(resp1)

figure; plot(resp1(1:5,:)')






end