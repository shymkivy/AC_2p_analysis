function f_plot_comp_scatter(score, ens_list, params, ops)

marker_type = f_get_param(params, 'marker_type', 'o');
marker_size = f_get_param(params, 'marker_size', 5);
mean_marker_type = f_get_param(params, 'mean_marker_type', '*');
mean_marker_size = f_get_param(params, 'mean_marker_size', 10);
plot_mean = f_get_param(params, 'plot_mean', 1);

[num_var, num_comp] = size(score);

if isempty(ens_list)
    ens_list = zeros(num_var,1);
end

ens_types = unique(ens_list);
ens_types(ens_types == 0) = [];

num_groups = numel(ens_types);
lg_ax = cell(num_groups,1);
lg_tit = cell(num_groups,1);
num_plt = 1;
figure; hold on;
if num_comp < 3
    gr_list = find(ens_list == 0);
    if ~isempty(gr_list)
        lg_ax{num_plt} = plot(score(gr_list,1), score(gr_list,end), marker_type, 'MarkerSize', marker_size, 'Color', 'k');
        lg_tit{num_plt} = 'core cells';
        num_plt = num_plt + 1;
    end
    if plot_mean
        gr_mean = mean(score(gr_list,:));
        lg_ax{num_groups+num_plt} = plot(gr_mean(1),gr_mean(2), mean_marker_type, 'MarkerSize', mean_marker_size, 'Color', 'k');
        lg_tit{num_groups+num_plt} = 'centeres';
    end
    for n_gr = 1:num_groups
        gr_list = find(ens_list == n_gr);
        lg_ax{num_plt+n_gr-1} = plot(score(gr_list,1),score(gr_list,2), marker_type, 'MarkerSize', marker_size, 'Color', ops.colors_list{n_gr});
        lg_tit{num_plt+n_gr-1} = sprintf('%d',n_gr);
        if plot_mean
            gr_mean = mean(score(gr_list,:));
            plot(gr_mean(1),gr_mean(2), mean_marker_type, 'MarkerSize', mean_marker_size, 'Color', ops.colors_list{n_gr});
        end
    end
else
    gr_list = find(ens_list == 0);
    if ~isempty(gr_list)
        lg_ax{num_plt} = plot3(score(gr_list,1),score(gr_list,2),score(gr_list,3), marker_type, 'MarkerSize', marker_size, 'Color', 'k');
        lg_tit{num_plt} = 'core cells';
        num_plt = num_plt + 1;
    end
    if plot_mean
        gr_mean = mean(score(gr_list,:));
        lg_ax{num_groups+num_plt} = plot3(gr_mean(1),gr_mean(2),gr_mean(3), mean_marker_type, 'MarkerSize', mean_marker_size, 'Color', 'k');
        lg_tit{num_groups+num_plt} = 'centeres';
    end
    for n_gr = 1:num_groups
        gr_list = find(ens_list == n_gr);
        lg_ax{num_plt+n_gr-1} = plot3(score(gr_list,1),score(gr_list,2),score(gr_list,3), marker_type, 'MarkerSize', marker_size, 'Color', ops.colors_list{n_gr});
        lg_tit{num_plt+n_gr-1} = sprintf('%d',n_gr);
        if plot_mean
            gr_mean = mean(score(gr_list,:));
            plot3(gr_mean(1),gr_mean(2),gr_mean(3), mean_marker_type, 'MarkerSize', mean_marker_size, 'Color', ops.colors_list{n_gr});
        end
    end
end
legend([lg_ax{:}], lg_tit)



end