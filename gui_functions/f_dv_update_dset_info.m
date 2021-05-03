function f_dv_update_dset_info(app)

app.mplSpinner.Value = min([app.mplSpinner.Value, app.ddata.num_planes]);

f_dv_update_A(app);
f_dv_update_trace(app);

end