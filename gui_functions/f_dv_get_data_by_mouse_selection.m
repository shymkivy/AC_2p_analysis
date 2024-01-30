function [data, title_tag1] = f_dv_get_data_by_mouse_selection(app)

button_text = app.SelectdatagroupDropDown.Value;

if strcmpi(button_text, 'plane')
    idx_data = strcmpi(app.data.dset_name_full, app.ddata.dset_name_full);
    title_tag = ['plane ' app.ddata.dset_name_full{1}];
elseif strcmpi(button_text, 'dataset')
    idx_data = strcmpi(app.data.dset_name_full, app.ddata.dset_name_full);
    title_tag = ['dset ' app.ddata.dset_name_full{1}];
elseif strcmpi(button_text, 'mouse region')
    idx_data = logical(strcmpi(app.data.mouse_id, app.ddata.mouse_id).*strcmpi(app.data.area, app.ddata.area));
    title_tag = ['mouse region  ' app.ddata.dset_name_full{1}];
elseif strcmpi(button_text, 'mouse')
    idx_data = strcmpi(app.data.mouse_id, app.ddata.mouse_id);
    title_tag = ['mouse ' app.ddata.mouse_id{1}];
elseif strcmpi(button_text, 'region')
    idx_data = strcmpi(app.data.area, app.ddata.area);
    title_tag = ['region ' app.ddata.area{1}];
elseif strcmpi(button_text, 'all')
    idx_data = true(numel(app.data.mouse_id),1);
    title_tag = 'all data';
end

data = app.data(idx_data,:);

title_tag1 = sprintf('%s; %s', data.paradigm{1}, title_tag);

end