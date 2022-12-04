function reg_cell_labels = f_dv_get_area_label(app, data1)

%num_cells = sum(cat(1,data1.stats{:}).num_cells);
num_cells = data1.num_cells;
reg_all = app.ops.regions_to_analyze;

% get area labels
if app.UseregdatalabelsCheckBox.Value
    if ~isempty(data1.registered_data{1})
        reg_cell_labels = data1.registered_data{1}.reg_labels;
    else
        reg_cell_labels = zeros(num_cells,1);
    end
else
    reg_idx = find(strcmpi(reg_all, data1.area));
    reg_cell_labels = ones(num_cells,1)*reg_idx;
end

end