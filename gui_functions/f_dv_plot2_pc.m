function f_dv_plot2_pc(data_in, tn_all, pcs, title_tag, color_cell, plot_idx)

if ~exist('pcs', 'var') || isempty(pcs)
    pcs = [1 2];
end

num_tn = numel(tn_all);

if ~exist('plot_idx', 'var') || isempty(plot_idx)
    plot_idx = {1:num_tn};
end

num_pl = numel(plot_idx);

figure; hold on
for n_pl = 1:num_pl
    plot(data_in(plot_idx{n_pl},pcs(1)), data_in(plot_idx{n_pl},pcs(2)), 'o-k', 'Linewidth', 1)
end
for n_tn = 1:num_tn
    %plot3([0 lr_data3d(n_tn, 1)], [0 lr_data3d(n_tn, 2)], [0 lr_data3d(n_tn, 3)], '-', 'color', app.ops.context_types_all_colors2{tn_all(n_tn)}, 'LineWidth', 2)
    plot(data_in(n_tn, pcs(1)), data_in(n_tn, pcs(2)), '.', 'color', color_cell{tn_all(n_tn)}, 'LineWidth', 2, 'MarkerSize', 20)
end
title(sprintf('2d %s', title_tag), 'interpreter', 'none');
xlabel('PC 1');
ylabel('PC 2');

end