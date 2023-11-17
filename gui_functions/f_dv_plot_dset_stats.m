function f_dv_plot_dset_stats(app)

data = app.data;

num_cells_all = cat(1,data.num_cells_pl{:});

fprintf('%d datasets; mean %.2f std %.2f sem %.2f cells per dset\n', numel(num_cells_all), mean(num_cells_all), std(num_cells_all), std(num_cells_all)/sqrt(numel(num_cells_all)-1))

end 
