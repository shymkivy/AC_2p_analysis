function f_dv_compute_stats(app)
% 
params = f_dv_gather_params(app);
params.ddata = app.ddata;
for n_pl = 1:params.ddata.num_planes
    params.n_pl = n_pl;
    params.cdata = app.cdata{n_pl};
    app.data(app.current_data_idx,:).stats{n_pl} = f_dv_compute_stats_core(app, params);
end
app.ddata = app.data(app.current_data_idx,:);

f_dv_update_dset_info(app);
f_dv_initialize_contours(app);
f_dv_set_contorus(app);

end