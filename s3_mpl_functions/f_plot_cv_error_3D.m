function f_plot_cv_error_3D(data, x_var, y_var, z_var)


data_x = [data.(x_var)];
data_y = [data.(y_var)];
data_z = [data.(z_var)];

numel_x = numel(unique(data_x));
numel_y = numel(unique(data_y));

grid_x = reshape(data_x(:), [numel_x, numel_y]);
grid_y = reshape(data_y(:), [numel_x, numel_y]);
grid_z = reshape(data_z(:), [numel_x, numel_y]);

figure; surf(grid_x, grid_y, grid_z)
title(z_var, 'interpreter', 'none')
xlabel(x_var, 'interpreter', 'none')
ylabel(y_var, 'interpreter', 'none')
end