function f_dv_plot2_pc_bytn(data_in, tn_all, tn_idx, pcs, title_tag, color_cell, render_painters)

if ~exist('pcs', 'var') || isempty(pcs)
    pcs = [1 2];
end

if ~exist('render_painters', 'var') || isempty(render_painters)
    render_painters = 0;
end

num_tn = numel(tn_all);

if render_painters
    figure(render='painters'); 
else
    figure();
end
hold on;
plot(data_in(:,pcs(1)), data_in(:,pcs(2)), 'ok', 'LineWidth', 2, 'MarkerSize', 6)

for n_tn = 1:num_tn
    idx1 = tn_idx == n_tn;
    %plot3([0 lr_data3d(n_tn, 1)], [0 lr_data3d(n_tn, 2)], [0 lr_data3d(n_tn, 3)], '-', 'color', app.ops.context_types_all_colors2{tn_all(n_tn)}, 'LineWidth', 2)
    plot(data_in(idx1, pcs(1)), data_in(idx1, pcs(2)), '.', 'color', color_cell{tn_all(n_tn)}, 'MarkerSize', 18)
end
title(sprintf('2d %s', title_tag), 'interpreter', 'none');
xlabel('PC 1');
ylabel('PC 2');

end