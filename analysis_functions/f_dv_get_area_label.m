function reg_cell_labels = f_dv_get_area_label(data, params, ops)

%num_cells = sum(cat(1,data1.stats{:}).num_cells);
num_cells = data.num_cells;
reg_all = ops.regions_to_analyze;

% get area labels
if params.use_reg_data_labels
    if ~isempty(data.registered_data{1})
        reg_cell_labels = data.registered_data{1}.reg_labels;
    else
        reg_cell_labels = zeros(num_cells,1);
    end
else
    reg_idx = find(strcmpi(reg_all, data.area));
    reg_cell_labels = ones(num_cells,1)*reg_idx;
end

end