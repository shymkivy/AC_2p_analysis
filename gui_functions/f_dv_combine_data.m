function [idx_stat, gr_tag] = f_dv_combine_data(app, mouse_id, dset_id)

gr_tag = app.statsbetweenDropDown.Value;
num_dsets = numel(mouse_id);

if strcmpi(gr_tag, 'Combined')
    idx_stat = ones(num_dsets,1);
elseif strcmpi(gr_tag, 'Mouse')
    mouse_idall = unique(mouse_id);
    idx_stat = zeros(num_dsets,1);
    for n_gr = 1:numel(mouse_idall)
        idx_stat(strcmpi(mouse_idall{n_gr}, mouse_id)) = n_gr;
    end
elseif strcmpi(gr_tag, 'Dset')
    idx_stat = dset_id;
elseif strcmpi(gr_tag, 'Subdset')
    idx_stat = (1:num_dsets)';
end

end