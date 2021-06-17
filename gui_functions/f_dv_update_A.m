function f_dv_update_A(app)

n_pl = app.mplSpinner.Value;

A = app.ddata.OA_data{n_pl}.est.A(:,app.cdata.accepted_cells);
A = A/max(A(:));

A_flat = reshape(full(sum(A,2)), 256, 256);

app.gui_plots.A_image.CData = A_flat;



end