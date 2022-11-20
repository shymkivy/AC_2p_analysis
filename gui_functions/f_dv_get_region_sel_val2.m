function [region_num, reg_tag] = f_dv_get_region_sel_val2(app)

reg_all = app.ops.regions_to_analyze;
reg_tag = app.regiontoplotDropDown.Value;
if strcmpi(reg_tag, 'All')
    region_num = (1:numel(reg_all))';
elseif strcmpi(reg_tag, 'All comb')
    region_num = 1:numel(reg_all);
elseif strcmpi(reg_tag, 'A1')
    region_num = 1;%find(strcmpi(reg_all, 'A1'));
elseif strcmpi(reg_tag, 'A2')
    region_num = 2;%find(strcmpi(reg_all, 'A2'));
elseif strcmpi(reg_tag, 'AAF')
    region_num = 3;%find(strcmpi(reg_all, 'AAF'));
elseif strcmpi(reg_tag, 'UF')
    region_num = 4;%find(strcmpi(reg_all, 'UF'));
elseif strcmpi(reg_tag, 'Primary vs secondary')
    region_num = [1, 3; 2, 4];
end

end