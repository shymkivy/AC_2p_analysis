function cell_plane_lut = f_dv_opto_make_cell_plane_lut(data)

num_cells = data.num_cells;

cell_plane_lut = [(1:num_cells)', zeros(num_cells,2)];
start1 = 1;
for n_pl = 1:data.num_planes
    num_cells2 = data.num_cells_pl{n_pl};
    end1 = start1 + num_cells2 - 1;
    cell_plane_lut(start1:end1,2) = n_pl;
    cell_plane_lut(start1:end1,3) = 1:num_cells2;
    start1 = end1 + 1;
end

end