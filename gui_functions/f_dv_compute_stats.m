function f_dv_compute_stats(app)
% 
params = f_dv_gather_params(app);

params.cdata = app.cdata;
params.ddata = app.ddata;
app.data(app.current_data_idx,:).stats{params.n_pl} = f_dv_compute_stats_core(app, params);

f_dv_update_dset_info(app);
f_dv_initialize_contours(app);
f_dv_set_contorus(app);

end