function f_dv_update_dset_info(app)

n_pl = app.mplSpinner.Value;

app.CellSpinner.Limits = [1, app.ddata.num_cells_pl{n_pl}];

f_dv_update_A(app);
f_dv_update_trace(app);

end