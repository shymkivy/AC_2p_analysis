function f_dv_stats_within(app)

params = f_dv_gather_params(app);
params.ddata = app.ddata;

fprintf('Computing stats _/%d planes: ', params.ddata.num_planes);
for n_pl = 1:params.ddata.num_planes
    fprintf('%d..',n_pl);
    params.n_pl = n_pl;
    params.cdata = app.cdata{n_pl};
    app.data(app.current_data_idx,:).stats_within{n_pl} = f_dv_stats_within_core(app, params);
end
app.ddata = app.data(app.current_data_idx,:);
fprintf('\nDone\n');

f_dv_update_dset_info(app);
f_dv_initialize_contours(app);
f_dv_set_contorus(app);
f_dv_update_params(app);

end