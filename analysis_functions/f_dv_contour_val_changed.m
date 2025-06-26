function f_dv_contour_val_changed(app)

min_val = app.ContourMinEditField.Value;
max_val = app.ContourMaxEditField.Value;

if min_val < app.gui_ops.contour_params.c_abs_lim(1)
    app.ContourMinEditField.Value = app.gui_ops.contour_params.c_abs_lim(1);
elseif min_val > app.ContourMaxEditField.Value
    app.ContourMinEditField.Value = max_val;
end

if max_val > app.gui_ops.contour_params.c_abs_lim(2)
    app.ContourMaxEditField.Value = app.gui_ops.contour_params.c_abs_lim(2);
elseif max_val < app.ContourMinEditField.Value
    app.ContourMaxEditField.Value = min_val;
end

app.gui_ops.contour_params.c_lim(1) = app.ContourMinEditField.Value;
app.gui_ops.contour_params.c_lim(2) = app.ContourMaxEditField.Value;
f_dv_update_contour_lims(app);
f_dv_set_contorus(app);

end