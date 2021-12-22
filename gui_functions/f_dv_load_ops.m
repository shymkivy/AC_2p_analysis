function f_dv_load_ops(app)

app.working_ops.trial_window = app.ops.trial_window;

% if app.ManualtimewinCheckBox.Value
%     app.working_ops.trial_window = [app.BaselineEditField.Value app.RespEditField.Value];
% else
%     app.working_ops.trial_window = app.ops.trial_window;
% end

end