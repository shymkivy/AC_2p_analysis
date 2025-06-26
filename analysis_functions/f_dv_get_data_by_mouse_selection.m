function [data, title_tag1] = f_dv_get_data_by_mouse_selection(data, params)

button_text = params.data_selection;

ddata = data(params.current_dset_idx,:);

if strcmpi(button_text, 'plane')
    %idx_data = strcmpi(app.data.dset_name_full, app.ddata.dset_name_full);
    idx_data = params.current_dset_idx;
    title_tag = ['plane ' ddata.dset_name_full{1}];
elseif strcmpi(button_text, 'dataset')
    idx_data = params.current_dset_idx;
    title_tag = ['dset ' ddata.dset_name_full{1}];
elseif strcmpi(button_text, 'mouse region')
    idx_data = logical(strcmpi(data.mouse_id, ddata.mouse_id).*strcmpi(data.area, ddata.area));
    title_tag = ['mouse region  ' ddata.dset_name_full{1}];
elseif strcmpi(button_text, 'mouse')
    idx_data = strcmpi(data.mouse_id, ddata.mouse_id);
    title_tag = ['mouse ' ddata.mouse_id{1}];
elseif strcmpi(button_text, 'region')
    idx_data = strcmpi(data.area, ddata.area);
    title_tag = ['region ' ddata.area{1}];
elseif strcmpi(button_text, 'all')
    idx_data = true(numel(data.mouse_id),1);
    title_tag = 'all data';
end

data = data(idx_data,:);

title_tag1 = sprintf('%s; %s', data.paradigm{1}, title_tag);

end