function f_plot_comp_scatter(score, ens_list, ops)

[num_var, num_comp] = size(score);

marker_type = 'o';
marker_size = 5;
mean_marker_type = '*';
mean_marker_size = 10;

if isempty(ens_list)
    ens_list = zeros(num_var,1);
end

ens_types = unique(ens_list);
ens_types(ens_types == 0) = [];

num_groups = numel(ens_types);
lg_ax = cell(num_groups+1,1);
lg_tit = cell(num_groups+1,1);
figure; hold on;
if num_comp < 3
    gr_list = find(ens_list == 0);
    lg_ax{1} = plot(score(gr_list,1), score(gr_list,end), marker_type, 'MarkerSize', marker_size, 'Color', 'k');
    lg_tit{1} = 'core cells';
    gr_mean = mean(score(gr_list,:));
    lg_ax{num_groups+2} = plot(gr_mean(1),gr_mean(2), mean_marker_type, 'MarkerSize', mean_marker_size, 'Color', 'k');
    lg_tit{num_groups+2} = 'centeres';
    for n_gr = 1:num_groups
        gr_list = find(ens_list == n_gr);
        lg_ax{n_gr+1} = plot(score(gr_list,1),score(gr_list,2), marker_type, 'MarkerSize', marker_size, 'Color', ops.colors_list{n_gr});
        lg_tit{n_gr+1} = sprintf('%d',n_gr);
        gr_mean = mean(score(gr_list,:));
        plot(gr_mean(1),gr_mean(2), mean_marker_type, 'MarkerSize', mean_marker_size, 'Color', ops.colors_list{n_gr});
    end
else
    gr_list = find(ens_list == 0);
    lg_ax{1} = plot3(score(gr_list,1),score(gr_list,2),score(gr_list,3), marker_type, 'MarkerSize', marker_size, 'Color', 'k');
    lg_tit{1} = 'core cells';
    gr_mean = mean(score(gr_list,:));
    lg_ax{num_groups+2} = plot3(gr_mean(1),gr_mean(2),gr_mean(3), mean_marker_type, 'MarkerSize', mean_marker_size, 'Color', 'k');
    lg_tit{num_groups+2} = 'centeres';
    for n_gr = 1:num_groups
        gr_list = find(ens_list == n_gr);
        lg_ax{n_gr+1} = plot3(score(gr_list,1),score(gr_list,2),score(gr_list,3), marker_type, 'MarkerSize', marker_size, 'Color', ops.colors_list{n_gr});
        lg_tit{n_gr+1} = sprintf('%d',n_gr);
        gr_mean = mean(score(gr_list,:));
        plot3(gr_mean(1),gr_mean(2),gr_mean(3), mean_marker_type, 'MarkerSize', mean_marker_size, 'Color', ops.colors_list{n_gr});
    end
end
legend([lg_ax{:}], lg_tit)



end