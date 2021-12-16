function f_dv_update_current_cdata(app)

params = f_dv_gather_params(app);

cdata_all = cell(app.ddata.num_planes,1);
for n_pl2 = 1:app.ddata.num_planes
    params.n_pl = n_pl2;
    cdata_all{n_pl2} = f_dv_compute_cdata(app, params);
end
app.cdata = cdata_all;

end