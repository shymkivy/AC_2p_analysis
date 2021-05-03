function f_dv_initialize(app)

app.A_image = imagesc(app.UIAxes, 0);
axis(app.UIAxes, 'equal');
axis(app.UIAxes, 'tight')

app.trace_plot = plot(app.UIAxes2, 0,0);
%axis(app.UIAxes2, 'tight');

end