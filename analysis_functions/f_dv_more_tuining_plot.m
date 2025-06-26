function f_dv_more_tuining_plot(app)

ops = app.ops;
params = f_dv_gather_params(app);

[data, title_tag] = f_dv_get_data_by_mouse_selection(app.data, params);


stats1 = data.stats{1};


resp1 = stats1.peak_vals(:,1:10);



figure; imagesc(resp1)

figure; plot(resp1(1:5,:)')






end