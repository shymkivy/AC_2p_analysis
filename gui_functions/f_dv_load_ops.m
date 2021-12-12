function f_dv_load_ops(app)

if app.ManualtimewinCheckBox.Value
    app.working_ops.trial_window = [app.BaselineEditField.Value, app.baserespsEditField.Value];
else
    app.working_ops.trial_window = app.ops.trial_window;
end

end