function [region_num, reg_tag, leg_list, reg_col] = f_dv_get_region_sel_val(params, ops)

reg_all = ops.regions_to_analyze;
reg_tag = params.region;
leg_list = {reg_tag};
reg_col = {'k'};
if strcmpi(reg_tag, 'All')
    region_num = (1:numel(reg_all))';
    leg_list = {'A1', 'A2', 'AAF', 'UF'};
    reg_col = ops.cond_colors;
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
    leg_list = {'Primary', 'Secondary'};
    reg_col = ops.cond_colors(1:2);
end

end