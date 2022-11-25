function color1 = f_dv_get_leg_color(app, leg1)

idx1 = strcmpi({app.ops.color_labels{:,1}}, leg1);
color1 = app.ops.color_labels{idx1,2};

end