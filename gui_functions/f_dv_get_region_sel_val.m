function [region_num, reg_tag] = f_dv_get_region_sel_val(app)

reg_all = app.ops.regions_to_analyze;
reg_tag = app.regiontoplotDropDown.Value;
if strcmpi(reg_tag, 'all')
    region_num = 1:numel(reg_all);
elseif strcmpi(reg_tag, 'A1')
    region_num = find(strcmpi(reg_all, 'A1'));
elseif strcmpi(reg_tag, 'A2')
    region_num = find(strcmpi(reg_all, 'A2'));
elseif strcmpi(reg_tag, 'AAF')
    region_num = find(strcmpi(reg_all, 'AAF'));
elseif strcmpi(reg_tag, 'UF')
    region_num = find(strcmpi(reg_all, 'UF'));
elseif strcmpi(reg_tag, 'All comb')
    region_num = (1:numel(reg_all))';
end

end