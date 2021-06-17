function [data, title_tag] = f_dv_get_data_by_mouse_selection(app)

button_text = app.SelectdatagroupButtonGroup.SelectedObject.Text;

if strcmpi(button_text, 'dset')
    idx_data = strcmpi(app.data.experiment, app.ddata.experiment);
    title_tag = ['dset ' app.ddata.experiment{1}];
elseif strcmpi(button_text, 'mouse')
    idx_data = strcmpi(app.data.mouse_tag, app.ddata.mouse_tag);
    title_tag = ['mouse ' app.ddata.mouse_tag{1}];
elseif strcmpi(button_text, 'region')
    idx_data = strcmpi(app.data.area, app.ddata.area);
    title_tag = ['region ' app.ddata.area{1}];
elseif strcmpi(button_text, 'all')
    idx_data = true(numel(app.data.experiment),1);
    title_tag = ['all data'];
end

data = app.data(idx_data,:);

end