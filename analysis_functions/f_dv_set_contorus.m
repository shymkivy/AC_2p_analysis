function f_dv_set_contorus(app)

params = f_dv_gather_params(app);
n_pl = params.n_pl;
ops = app.ops;

accepted_cells = app.cdata{n_pl}.accepted_cells;
num_cells = app.cdata{n_pl}.num_cells;

use_color_map = 1;

stats1 = app.ddata.stats{n_pl};

app.gui_ops.contour_params.visible_set = 0;
app.gui_ops.contour_params.contour_mag = zeros(num_cells,1);
app.gui_ops.contour_params.c_lim = [0 0];
app.gui_ops.contour_params.c_abs_lim = [0, 0];

contour_val = app.ContoursDropDown.Value;

if strcmpi(contour_val, 'SNR')
    app.gui_ops.contour_params.visible_set = 1;
    SNR_list = app.ddata.OA_data{n_pl}.proc.SNR2_vals(accepted_cells);
    app.gui_ops.contour_params.contour_mag = SNR_list; % factor to resize magnitudes to fir color
    app.gui_ops.contour_params.c_abs_lim = [min(SNR_list) max(SNR_list)];
    if isfield(app.gui_ops, 'SNR_lim')
        app.gui_ops.contour_params.c_lim = app.gui_ops.SNR_lim;
    else
        app.gui_ops.contour_params.c_lim = app.gui_ops.contour_params.c_abs_lim;
    end
elseif strcmpi(contour_val, 'Skewness')
    app.gui_ops.contour_params.visible_set = 1;
    skew_list = app.ddata.OA_data{n_pl}.proc.skewness(accepted_cells);
    app.gui_ops.contour_params.contour_mag = skew_list; % factor to resize magnitudes to fir color
    app.gui_ops.contour_params.c_abs_lim = [min(skew_list) max(skew_list)];
    if isfield(app.gui_ops, 'skew_lim')
        app.gui_ops.contour_params.c_lim = app.gui_ops.skew_lim;
    else
        app.gui_ops.contour_params.c_lim = app.gui_ops.contour_params.c_abs_lim;
    end
elseif ~isempty(stats1)
    if strcmpi(contour_val, 'Tuning type')
        tn_all = f_dv_get_trial_number(params);
        
        [resp_cells, ~, resp_vals] = f_dv_get_resp_vals_cells(stats1, tn_all, params);

        resp_vals(~resp_cells) = 0;
        [max_val, max_idx] = max(resp_vals, [], 2);
        app.gui_ops.contour_params.visible_set = logical(max_val);
        app.gui_ops.contour_params.contour_mag = max_idx;
        app.gui_ops.contour_params.c_abs_lim = [min(max_idx) max(max_idx)];
        app.gui_ops.contour_params.c_lim = [min(max_idx) max(max_idx)];

        use_color_map = 0;
    elseif strcmpi(contour_val, 'Tuning magnitude')
        tn_all = f_dv_get_trial_number(params);
        
        [resp_cells, ~, resp_vals] = f_dv_get_resp_vals_cells(stats1, tn_all, params);
        
        if app.ConverttoZCheckBox.Value
            st_mean_mean = stats1.stat_trials_mean_mean;
            st_mean_sem = stats1.stat_trials_mean_sem;
            resp_vals = (resp_vals - st_mean_mean)./st_mean_sem;
        end
        resp_vals(~resp_cells) = 0;
        resp_vals2 = max(resp_vals,[],2);
        app.gui_ops.contour_params.visible_set = resp_vals2>0;

        app.gui_ops.contour_params.contour_mag = resp_vals2; % factor to resize magnitudes to fir color
        app.gui_ops.contour_params.c_abs_lim = [min(resp_vals2) max(resp_vals2)];

        if isfield(app.gui_ops, 'tuning_lim')
            app.gui_ops.contour_params.c_lim = app.gui_ops.tuning_lim;
        else
            app.gui_ops.contour_params.c_lim = app.gui_ops.contour_params.c_abs_lim;
        end
    end    
end
visible_set = app.gui_ops.contour_params.visible_set;

if use_color_map
    c_lim = app.gui_ops.contour_params.c_lim;
    c_lim(2) = min(c_lim(2), 100);
    app.ContourMinEditField.Value = c_lim(1);
    app.ContourMaxEditField.Value = c_lim(2);

    % crop if goes beyond lim
    contour_mag = max(app.gui_ops.contour_params.contour_mag, c_lim(1));
    contour_mag = min(contour_mag, c_lim(2));

    % create color map
    contour_resolution = 50;
    color_map = jet(round(contour_resolution*c_lim(2)) - round(contour_resolution*c_lim(1))+1); 
    color_index = round(contour_resolution*c_lim(1)):round(contour_resolution*c_lim(2));
    imagesc(app.UIAxesColorPallet, color_index/contour_resolution, 1, reshape(color_map,1,[],3));
    axis(app.UIAxesColorPallet, 'tight');
else
    contour_mag = app.gui_ops.contour_params.contour_mag;
    imagesc(app.UIAxesColorPallet, reshape([ops.context_types_all_colors2{1:10}]', 1, 10, 3));
end

% update colors
if numel(visible_set)>1
    for n_cell = 1:num_cells
        app.gui_plots.contours_gobj(n_cell).Visible = visible_set(n_cell);
        if use_color_map
            app.gui_plots.contours_gobj(n_cell).Color = color_map(round(contour_resolution*contour_mag(n_cell)) == color_index,:); 
        else
            app.gui_plots.contours_gobj(n_cell).Color = ops.context_types_all_colors2{tn_all(contour_mag(n_cell))};
        end
    end
else
    for n_cell = 1:num_cells
        app.gui_plots.contours_gobj(n_cell).Visible = visible_set;
        if use_color_map
            app.gui_plots.contours_gobj(n_cell).Color = color_map(round(contour_resolution*contour_mag(n_cell)) == color_index,:); 
        else
            app.gui_plots.contours_gobj(n_cell).Color = ops.context_types_all_colors2{tn_all(contour_mag(n_cell))};
        end
    end
end


end