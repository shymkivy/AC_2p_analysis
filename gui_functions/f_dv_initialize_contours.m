function f_dv_initialize_contours(app)

for n_obj = 1:numel(app.gui_plots.contours_gobj)
    delete(app.gui_plots.contours_gobj(n_obj))
end

n_pl = app.mplSpinner.Value;

est1 = app.ddata.OA_data{n_pl}.est;

cdata = app.cdata{n_pl};
num_cells = cdata.num_cells;
cell_num_convert = find(cdata.accepted_cells);

app.gui_plots.contours_gobj = gobjects(num_cells,1);
hold(app.UIAxes, 'on');
for n_cell = 1:num_cells
    temp_contours = est1.contours{cell_num_convert(n_cell)};
    app.gui_plots.contours_gobj(n_cell) = plot(app.UIAxes, temp_contours(:,2), temp_contours(:,1), 'LineWidth', 1, 'Visible', 0);
end
hold(app.UIAxes, 'off');

end