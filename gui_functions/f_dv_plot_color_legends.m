function f_dv_plot_color_legends(app)

map1 = reshape(cat(1,app.ops.context_types_all_colors2{1:10}), [10, 1, 3]);

figure; imagesc(map1); axis equal tight
ylabel('freq')

end