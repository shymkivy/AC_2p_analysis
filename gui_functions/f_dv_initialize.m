function f_dv_initialize(app)

app.gui_plots.A_image = imagesc(app.UIAxes, 0);
app.gui_plots.A_image.ButtonDownFcn = @(~,~) f_dv_button_down(app, app.gui_plots.A_image);
app.gui_plots.plot_current_contour = [];
axis(app.UIAxes, 'equal');
axis(app.UIAxes, 'tight')

hold(app.UIAxes2, 'on');
app.gui_plots.plot_raw = plot(app.UIAxes2, 0,0, 'color', [0.8500, 0.3250, 0.0980]);
app.gui_plots.plot_C = plot(app.UIAxes2, 0,0, 'color', [0, 0.4470, 0.7410]);
app.gui_plots.plot_spikes = plot(app.UIAxes2, 0,0, 'color', [0.4940, 0.1840, 0.5560]);
app.gui_plots.plot_stim_times = plot(app.UIAxes2, 0,0, 'color', [0, 0, 0]);
%axis(app.UIAxes2, 'tight');

app.gui_plots.freq_resp_fig = [];
app.gui_plots.ctx_resp_fig = [];
app.gui_plots.select_resp_fig = [];
app.gui_plots.contours_gobj = [];

app.StatsourceDropDown.Items = {'Freqs', 'All'};
app.StatmethodDropDown.Items = {'Sample', 'Pop percentile'};

end