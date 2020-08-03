function f_plot_dset_deets(plot_y_cell, ops)
figure; hold on;
xpts = numel(ops.regions_to_analyze);

lg1 = cell(xpts,1);
pl_n = cell(xpts,1);
for n_cond = 1:xpts
    num_x = numel(plot_y_cell{n_cond});
    pl_n_temp = plot(linspace(n_cond-0.1,n_cond+0.1,num_x), plot_y_cell{n_cond}, 'o', 'color', ops.cond_colors{n_cond});
    pl_n{n_cond} = pl_n_temp(1);
    mean_cond = mean(plot_y_cell{n_cond});
    line([n_cond-0.25 n_cond+0.25], [mean_cond mean_cond], 'Color', 'k', 'LineWidth', 2)
    lg1{n_cond} = ops.regions_to_analyze{n_cond};
end
legend([pl_n{:}], ops.regions_to_analyze);
pl = gca;
pl.XTickLabel = [];
pl.XTick = 1:xpts;
end