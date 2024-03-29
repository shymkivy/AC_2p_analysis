function f_dv_plot3_t_pc(top_comp, tn_all, title_tag, colors_tn, add_shadow, shadow_axis_locs, render_painters, add_grid, reverse_xyz)

if ~exist('add_shadow', 'var') || isempty(add_shadow)
    add_shadow = 1;
end

if ~exist('shadow_axis_locs', 'var') || isempty(shadow_axis_locs)
    shadow_axis_locs = [1 2 1];
end

shadow_alpha = 0.3;
shadow_line_width = 1;
%shadow_lim_offset = 1.5;
shadow_lim_offset = 0.5; %1.7, 1.5
non_shadow_lim_offset = 0.2; %1.4 1.3

[num_t, num_tn, ~] = size(top_comp);

if render_painters
    f1 = figure(render='painters'); 
else
    f1 = figure();
end
hold on; %axis padded
for n_tn = 1:num_tn
    tn1 = tn_all(n_tn);
    if sum(tn1 == [28, 29, 30])
        symb1 = '*';
    else
        symb1 = 'o';
    end
    plot3(top_comp(:,n_tn, 1), top_comp(:,n_tn, 2), top_comp(:,n_tn, 3), 'color', colors_tn{tn1}, 'LineWidth', 2);
    plot3(top_comp(1,n_tn, 1), top_comp(1,n_tn, 2), top_comp(1,n_tn, 3), symb1, 'color', colors_tn{tn1}, 'LineWidth', 2);
end
xlabel('PC 1');
ylabel('PC 2');
zlabel('PC 3');

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

    for n_tn = 1:num_tn
        tn1 = tn_all(n_tn);
%             if sum(tn1 == [28, 29, 30])
%                 symb1 = '*';
%             else
%                 symb1 = 'o';
%             end
        plot3(ones(num_t,1)*sh(1), top_comp(:,n_tn, 2), top_comp(:,n_tn, 3), 'color', [colors_tn{tn1} shadow_alpha], 'LineWidth', shadow_line_width);
        %pl1 = plot3(xyzlims(1), top_comp2(1,n_tn, trs1(2)), top_comp2(1,n_tn, trs1(3)), symb1, 'color', [colors_tn{tn1} shadow_alpha], 'LineWidth', shadow_line_width);
        plot3(top_comp(:,n_tn, 1), ones(num_t,1)*sh(2), top_comp(:,n_tn, 3), 'color',[colors_tn{tn1} shadow_alpha], 'LineWidth', shadow_line_width);
        %pl1 = plot3(top_comp2(1,n_tn, 1), xyzlims(2), top_comp2(1,n_tn, 3), symb1, 'color', [colors_tn{tn1} alpha1], 'LineWidth', shadow_line_width);
        plot3(top_comp(:,n_tn, 1), top_comp(:,n_tn, 2), ones(num_t,1)*sh(3), 'color', [colors_tn{tn1} shadow_alpha], 'LineWidth', shadow_line_width);
        %pl1 = plot3(top_comp2(1,n_tn, 1), top_comp2(1,n_tn, 2), xyzlims(3), symb1, 'color', [colors_tn{tn1} alpha1], 'LineWidth', shadow_line_width);
    end

    f1.Children(1).XLim = xyz_lh(1,:);
    f1.Children(1).YLim = xyz_lh(2,:);
    f1.Children(1).ZLim = xyz_lh(3,:);

end
title(title_tag, 'Interpreter','none');
grid on;
  
end