function f_plot_dset_deets(plot_y_cell, ops, sp_pt)
if ~exist('sp_pt', 'var')
    figure;
else
    subplot(sp_pt)
end
hold on;
xpts = numel(ops.regions_to_analyze);
plot_shuff = size(plot_y_cell,2)-1;
lg1 = cell(xpts,1);
pl_n = cell(xpts,1);
for n_cond = 1:xpts
    num_x = numel(plot_y_cell{n_cond,1});
    pl_n_temp = plot(linspace(n_cond-0.1,n_cond+0.1,num_x), plot_y_cell{n_cond,1}, 'o', 'color', ops.cond_colors{n_cond});
    pl_n{n_cond} = pl_n_temp(1);
    mean_cond = mean(plot_y_cell{n_cond,1});
    line([n_cond-0.25 n_cond+0.25], [mean_cond mean_cond], 'Color', 'k', 'LineWidth', 2);
    lg1{n_cond} = ops.regions_to_analyze{n_cond};
    if  plot_shuff
        plot(linspace(n_cond-0.1,n_cond+0.1,num_x), plot_y_cell{n_cond,2}, 'o', 'color', 'k');
    end
end
legend([pl_n{:}], ops.regions_to_analyze);
pl = gca;
pl.XTickLabel = [];
pl.XTick = 1:xpts;
end