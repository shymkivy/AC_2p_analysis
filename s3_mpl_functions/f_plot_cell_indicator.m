function f_plot_cell_indicator(raster, cell_list, ops)
hold on;
[num_cells, num_row] = size(raster);

cell_types = unique(cell_list);
cell_types(cell_types == 0) = [];

num_ct = numel(cell_types);

color_seq_tt = zeros(num_cells,1,3);
for n_ct = 1:num_ct
    cells1 = find(cell_list == cell_types(n_ct));
    num_list_cells1 = numel(cells1);
    for n_cell_ind = 1:num_list_cells1
        n_cell = cells1(n_cell_ind);
        color_seq_tt(n_cell,:,:) = ops.colors_list{n_ct};
    end
end

col_width = ceil(num_row/50);

imagesc(num_row+(1:col_width),1:num_cells,repmat(color_seq_tt,1,col_width,1));
axis tight

end