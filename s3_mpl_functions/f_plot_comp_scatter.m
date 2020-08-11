function f_plot_comp_scatter(score, ens_list, ops)

[~, num_comp] = size(score);

marker_size = 5;
marker_type = 'o';

ens_types = unique(ens_list);
ens_types(ens_types == 0) = [];

num_groups = numel(ens_types);
lg_ax = cell(num_groups,1);
lg_tit = cell(num_groups,1);
figure; hold on;
if num_comp < 3
    plot(score(:,1), score(:,end), marker_type, 'MarkerSize', marker_size, 'Color', 'k');
    for n_gr = 1:num_groups
        gr_list = find(ens_list == n_gr);
        lg_ax{n_gr} = plot(score(gr_list,1),score(gr_list,2), marker_type, 'MarkerSize', marker_size, 'Color', ops.colors_list{n_gr});
        lg_tit{n_gr} = sprintf('%d',n_gr);
    end
else
    plot3(score(:,1),score(:,2),score(:,3), marker_type, 'MarkerSize', marker_size, 'Color', 'k');
    for n_gr = 1:num_groups
        gr_list = find(ens_list == n_gr);
        lg_ax{n_gr} = plot3(score(gr_list,1),score(gr_list,2),score(gr_list,3), marker_type, 'MarkerSize', marker_size, 'Color', ops.colors_list{n_gr});
        lg_tit{n_gr} = sprintf('%d',n_gr);
    end
end
legend([lg_ax{:}], lg_tit)



end