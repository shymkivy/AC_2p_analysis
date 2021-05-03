function f_dv_random_button(app)

n_pl = app.mplSpinner.Value;
n_cell = app.CellSpinner.Value;

if n_pl>1
   num_cells_pl2 = app.ddata.num_cells_pl(1:(n_pl-1));
   n_cell_full = sum(cat(1,num_cells_pl2{:})) + n_cell;
else
    n_cell_full = n_cell;
end

1

end