function f_dv_update_trace(app)

n_pl = app.mplSpinner.Value;
n_cell = app.CellSpinner.Value;

app.trace_plot.XData = app.ddata.proc_data{1}.frame_data.frame_times_mpl{n_pl};
if strcmpi(app.ButtonGroup.SelectedObject.Text, 'Raw')
    app.trace_plot.YData = app.ddata.traces_raw{n_pl}(n_cell,:);
elseif strcmpi(app.ButtonGroup.SelectedObject.Text, 'Firing Rate')
    app.trace_plot.YData = app.ddata.firing_rate{n_pl}(n_cell,:);
elseif strcmpi(app.ButtonGroup.SelectedObject.Text, 'Smooth FR')
    app.trace_plot.YData = app.ddata.firing_rate_sm{n_pl}(n_cell,:);
end
end