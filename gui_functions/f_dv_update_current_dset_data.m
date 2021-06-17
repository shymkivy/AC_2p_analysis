function f_dv_update_current_dset_data(app)

params = f_dv_gather_params(app);
cdata = f_dv_compute_cdata(app, params);

app.cdata = cdata;

end