function color1 = f_dv_get_leg_color(ops, leg1)

idx1 = strcmpi({ops.color_labels{:,1}}, leg1);
color1 = ops.color_labels{idx1,2};

end