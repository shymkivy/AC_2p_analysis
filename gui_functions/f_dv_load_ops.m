function f_dv_load_ops(app)

if app.ManualtimewinCheckBox.Value
   try
        win1 = app.baserespsEditField.Value;
        k = strfind(app.baserespsEditField.Value,',');
        app.working_ops.trial_window = [str2double(win1(1:k-1)), str2double(win1(k+1:end))];
    catch
        app.working_ops.trial_window = app.ops.trial_window;
        fprintf('Cannot parse manual wime win input, using [%.2f, %.2f]\n', app.ops.trial_window(1), app.ops.trial_window(2));
    end
else
    app.working_ops.trial_window = app.ops.trial_window;
end

end