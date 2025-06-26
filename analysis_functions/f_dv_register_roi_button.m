function f_dv_register_roi_button(app)

disp('registering rois...')

ddata = app.ddata;

%[data, title_tag] = f_dv_get_data_by_mouse_selection(app);

idx1 = logical(strcmpi(ddata.mouse_id, app.data.mouse_id).*(ddata.FOV_num == app.data.FOV_num));
data2 = app.data(idx1,:);

app.data(idx1,:).register_roi = f_dv_register_roi_core(data2, app.regscoretreshEditField.Value);

disp('Done')

end