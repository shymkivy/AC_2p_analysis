function f_dv_initialize(app)

app.A_image = imagesc(app.UIAxes, 0);
app.A_image.ButtonDownFcn = @(~,~) f_dv_button_down(app, app.A_image);
axis(app.UIAxes, 'equal');
axis(app.UIAxes, 'tight')

hold(app.UIAxes2, 'on');
app.plot_raw = plot(app.UIAxes2, 0,0, 'color', [0.8500, 0.3250, 0.0980]);
app.plot_C = plot(app.UIAxes2, 0,0, 'color', [0, 0.4470, 0.7410]);
app.plot_spikes = plot(app.UIAxes2, 0,0, 'color', [0.4940, 0.1840, 0.5560]);
app.plot_stim_times = plot(app.UIAxes2, 0,0, 'color', [0, 0, 0]);
%axis(app.UIAxes2, 'tight');

end