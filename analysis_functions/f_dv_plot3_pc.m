function f_dv_plot3_pc(data_in, tn_all, pcs, title_tag, color_cell, plot_idx, add_shadow, shadow_axis_locs, render_painters, add_grid, reverse_xyz)

if ~exist('pcs', 'var') || isempty(pcs)
    pcs = [1 2 3];
end

if ~exist('add_shadow', 'var') || isempty(add_shadow)
    add_shadow = 1;
end

if ~exist('shadow_axis_locs', 'var') || isempty(shadow_axis_locs)
    shadow_axis_locs = [1 2 1];
end

if ~exist('render_painters', 'var') || isempty(render_painters)
    render_painters = 0;
end

if ~exist('reverse_xyz', 'var') || isempty(reverse_xyz)
    reverse_xyz = [0, 0, 0];
end

if ~exist('add_grid', 'var') || isempty(add_grid)
    add_grid = 1;
end

shadow_alpha = 0.2;
shadow_line_width = 1;
shadow_lim_offset = 0.5; %1.7, 1.5
non_shadow_lim_offset = 0.5; %1.4 1.3

num_tn = numel(tn_all);

if ~exist('plot_idx', 'var') || isempty(plot_idx)
    plot_idx = {1:num_tn};
end
num_pl = numel(plot_idx);

if render_painters
    f1 = figure(render='painters'); 
else
    f1 = figure();
end
hold on;
plot3(data_in(:,pcs(1)), data_in(:,pcs(2)), data_in(:,pcs(3)),'ok', 'LineWidth', 2, 'MarkerSize', 7);
for n_pl = 1:num_pl
    plot3(data_in(plot_idx{n_pl},pcs(1)), data_in(plot_idx{n_pl},pcs(2)), data_in(plot_idx{n_pl},pcs(3)), '-k', 'Linewidth', 1)
end

for n_tn = 1:num_tn
    %plot3([0 lr_data3d(n_tn, 1)], [0 lr_data3d(n_tn, 2)], [0 lr_data3d(n_tn, 3)], '-', 'color', app.ops.context_types_all_colors2{tn_all(n_tn)}, 'LineWidth', 2)
    plot3(data_in(n_tn, pcs(1)), data_in(n_tn, pcs(2)), data_in(n_tn, pcs(3)), '.', 'color', color_cell{tn_all(n_tn)}, 'MarkerSize', 20)
end

if add_shadow
    axis tight;
    sign_idx = [-1, 1];
    xyz_lh = [f1.Children(1).XLim; f1.Children(1).YLim; f1.Children(1).ZLim];
    xyz_mag = diff(xyz_lh,[],2);
    
    sh = zeros(3,1);
    for n1 = 1:3
        xyz_lh(n1, shadow_axis_locs(n1)) = xyz_lh(n1, shadow_axis_locs(n1)) + sign_idx(shadow_axis_locs(n1))*xyz_mag(n1)*shadow_lim_offset;
        xyz_lh(n1, 3-shadow_axis_locs(n1)) = xyz_lh(n1, 3-shadow_axis_locs(n1)) + sign_idx(3-shadow_axis_locs(n1))*xyz_mag(n1)*non_shadow_lim_offset;
        sh(n1) = xyz_lh(n1,shadow_axis_locs(n1));
    end

    for n_pl = 1:num_pl
        num_t = numel(plot_idx{n_pl});
        plot3(ones(num_t,1)*sh(1), data_in(plot_idx{n_pl}, pcs(2)), data_in(plot_idx{n_pl}, pcs(3)), 'color', [0 0 0 shadow_alpha], 'LineWidth', shadow_line_width);
        plot3(data_in(plot_idx{n_pl}, pcs(1)), ones(num_t,1)*sh(2), data_in(plot_idx{n_pl}, pcs(3)), 'color', [0 0 0 shadow_alpha], 'LineWidth', shadow_line_width);
        plot3(data_in(plot_idx{n_pl}, pcs(1)), data_in(plot_idx{n_pl}, pcs(2)), ones(num_t,1)*sh(3), 'color', [0 0 0 shadow_alpha], 'LineWidth', shadow_line_width);
    end

    f1.Children(1).XLim = xyz_lh(1,:);
    f1.Children(1).YLim = xyz_lh(2,:);
    f1.Children(1).ZLim = xyz_lh(3,:);

    for n_tn = 1:num_tn
        tn1 = tn_all(n_tn);
        sc1 = scatter3(sh(1), data_in(n_tn, pcs(2)), data_in(n_tn, pcs(3)), 314);
        sc1.Marker = '.';
        sc1.MarkerEdgeColor = color_cell{tn1};
        sc1.MarkerEdgeAlpha = shadow_alpha;
        sc1 = scatter3(data_in(n_tn, pcs(1)), sh(2), data_in(n_tn, pcs(3)), 314);
        sc1.Marker = '.';
        sc1.MarkerEdgeColor = color_cell{tn1};
        sc1.MarkerEdgeAlpha = shadow_alpha;
        sc1 = scatter3(data_in(n_tn, pcs(1)), data_in(n_tn, pcs(2)), sh(3), 314);
        sc1.Marker = '.';
        sc1.MarkerEdgeColor = color_cell{tn1};
        sc1.MarkerEdgeAlpha = shadow_alpha;
    end
end
title(sprintf('3d %s', title_tag), 'interpreter', 'none');
xlabel(sprintf('PC %d', pcs(1)));
ylabel(sprintf('PC %d', pcs(2)));
zlabel(sprintf('PC %d', pcs(3)));
if add_grid
    grid on;
end

if reverse_xyz(1)
    f1.Children(1).XDir = 'reverse';
end
if reverse_xyz(2)
    f1.Children(1).YDir = 'reverse';
end
if reverse_xyz(3)
    f1.Children(1).ZDir = 'reverse';
end

end