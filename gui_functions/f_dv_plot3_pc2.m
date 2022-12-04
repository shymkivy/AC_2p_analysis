function f_dv_plot3_pc2(data_in, tn_all, pcs, title_tag, color_cell, add_shadow, plot_idx)

if ~exist('pcs', 'var') || isempty(pcs)
    pcs = [1 2 3];
end

if ~exist('add_shadow', 'var') || isempty(add_shadow)
    add_shadow = 1;
end

shadow_alpha = 0.3;
shadow_line_width = 1;
shadow_lim_oddset = 1.5;
shadow_axis_locs = [1 2 1];

num_tn = numel(tn_all);

if ~exist('plot_idx', 'var') || isempty(plot_idx)
    plot_idx = {1:num_tn};
end
num_pl = numel(plot_idx);

f1 = figure; hold on
for n_pl = 1:num_pl
    plot3(data_in(plot_idx{n_pl},pcs(1)), data_in(plot_idx{n_pl},pcs(2)), data_in(plot_idx{n_pl},pcs(3)), 'o-k', 'Linewidth', 1)
end
for n_tn = 1:num_tn
    %plot3([0 lr_data3d(n_tn, 1)], [0 lr_data3d(n_tn, 2)], [0 lr_data3d(n_tn, 3)], '-', 'color', app.ops.context_types_all_colors2{tn_all(n_tn)}, 'LineWidth', 2)
    plot3(data_in(n_tn, pcs(1)), data_in(n_tn, pcs(2)), data_in(n_tn, pcs(3)), '.', 'color', color_cell{tn_all(n_tn)}, 'LineWidth', 2, 'MarkerSize', 20)
end
if add_shadow
    xyzlims = [f1.Children(1).XLim(shadow_axis_locs(1)), f1.Children(1).YLim(shadow_axis_locs(2)), f1.Children(1).ZLim(shadow_axis_locs(3))]*shadow_lim_oddset;
    
    for n_pl = 1:num_pl
        num_t = numel(plot_idx{n_pl});
        plot3(ones(num_t,1)*xyzlims(1), data_in(plot_idx{n_pl}, pcs(2)), data_in(plot_idx{n_pl}, pcs(3)), 'color', [0 0 0 shadow_alpha], 'LineWidth', shadow_line_width);
        plot3(data_in(plot_idx{n_pl}, pcs(1)), ones(num_t,1)*xyzlims(2), data_in(plot_idx{n_pl}, pcs(3)), 'color', [0 0 0 shadow_alpha], 'LineWidth', shadow_line_width);
        plot3(data_in(plot_idx{n_pl}, pcs(1)), data_in(plot_idx{n_pl}, pcs(2)), ones(num_t,1)*xyzlims(3), 'color', [0 0 0 shadow_alpha], 'LineWidth', shadow_line_width);
    end

    for n_tn = 1:num_tn
        tn1 = tn_all(n_tn);
        sc1 = scatter3(xyzlims(1), data_in(n_tn, pcs(2)), data_in(n_tn, pcs(3)), 314);
        sc1.Marker = '.';
        sc1.MarkerEdgeColor = color_cell{tn1};
        sc1.MarkerEdgeAlpha = shadow_alpha;
        sc1 = scatter3(data_in(n_tn, pcs(1)), xyzlims(2), data_in(n_tn, pcs(3)), 314);
        sc1.Marker = '.';
        sc1.MarkerEdgeColor = color_cell{tn1};
        sc1.MarkerEdgeAlpha = shadow_alpha;
        sc1 = scatter3(data_in(n_tn, pcs(1)), data_in(n_tn, pcs(2)), xyzlims(3), 314);
        sc1.Marker = '.';
        sc1.MarkerEdgeColor = color_cell{tn1};
        sc1.MarkerEdgeAlpha = shadow_alpha;
    end
end
title(sprintf('3d %s', title_tag), 'interpreter', 'none');
xlabel('PC 1');
ylabel('PC 2');
zlabel('PC 3');
grid on

end