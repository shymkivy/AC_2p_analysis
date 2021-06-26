function f_dv_set_contorus(app)

n_pl = app.mplSpinner.Value;

accepted_cells = app.cdata.accepted_cells;
num_cells = app.cdata.num_cells;

use_color_map = 1;

stats1 = app.ddata.stats{n_pl};

if app.ConverttoZCheckBox.Value
    pop_mean_val = stats1.pop_mean_val;
    pop_z_factor = stats1.pop_z_factor;
end

if strcmp(app.ContoursButtonGroup.SelectedObject.Text,'None')
    app.gui_ops.contour_params.visible_set = 0;
    app.gui_ops.contour_params.contour_mag = zeros(num_cells,1);
    app.gui_ops.contour_params.c_lim = [0 0];
    app.gui_ops.contour_params.c_abs_lim = [0, 0];
elseif strcmp(app.ContoursButtonGroup.SelectedObject.Text,'Tuning type')
    tn_all = f_dv_get_trial_number(app);
    tuning_freq = stats1.peak_val_all(:,tn_all);
    resp_cells = stats1.cell_is_resp(:,tn_all);
    tuning_freq(~resp_cells) = 0;
    [max_val, max_idx] = max(tuning_freq, [], 2);
    app.gui_ops.contour_params.visible_set = logical(max_val);
    app.gui_ops.contour_params.contour_mag = max_idx;
    app.gui_ops.contour_params.c_abs_lim = [min(max_idx) max(max_idx)];
    app.gui_ops.contour_params.c_lim = [min(max_idx) max(max_idx)];
    use_color_map = 0;
elseif strcmp(app.ContoursButtonGroup.SelectedObject.Text,'Tuning mag')
    tn_all = f_dv_get_trial_number(app);
    resp_cells = stats1.cell_is_resp;
    peak_vals = stats1.peak_val_all;
    if app.ConverttoZCheckBox.Value
        peak_vals = (peak_vals - pop_mean_val)./pop_z_factor;
    end
    peak_vals(~resp_cells) = 0;
    peak_vals2 = max(peak_vals(:,tn_all),[],2);
    app.gui_ops.contour_params.visible_set = peak_vals2>0;

    app.gui_ops.contour_params.contour_mag = peak_vals2; % factor to resize magnitudes to fir color
    app.gui_ops.contour_params.c_abs_lim = [min(peak_vals2) max(peak_vals2)];
    if isfield(app.gui_ops, 'tuning_lim')
        app.gui_ops.contour_params.c_lim = app.gui_ops.tuning_lim;
    else
        app.gui_ops.contour_params.c_lim = app.gui_ops.contour_params.c_abs_lim;
    end
elseif strcmp(app.ContoursButtonGroup.SelectedObject.Text,'SNR')
    app.gui_ops.contour_params.visible_set = 1;
    SNR_list = app.ddata.OA_data{n_pl}.proc.SNR2_vals(accepted_cells);
    app.gui_ops.contour_params.contour_mag = SNR_list; % factor to resize magnitudes to fir color
    app.gui_ops.contour_params.c_abs_lim = [min(SNR_list) max(SNR_list)];
    if isfield(app.gui_ops, 'SNR_lim')
        app.gui_ops.contour_params.c_lim = app.gui_ops.SNR_lim;
    else
        app.gui_ops.contour_params.c_lim = app.gui_ops.contour_params.c_abs_lim;
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
    imagesc(app.UIAxesColorPallet, []);
end

% update colors
if numel(visible_set)>1
    for n_cell = 1:num_cells
        app.gui_plots.contours_gobj(n_cell).Visible = visible_set(n_cell);
        if use_color_map
            app.gui_plots.contours_gobj(n_cell).Color = color_map(round(contour_resolution*contour_mag(n_cell)) == color_index,:); 
        else
            app.gui_plots.contours_gobj(n_cell).Color = app.ops.context_types_all_colors2{tn_all(contour_mag(n_cell))};
        end
    end
else
    for n_cell = 1:num_cells
        app.gui_plots.contours_gobj(n_cell).Visible = visible_set;
        if use_color_map
            app.gui_plots.contours_gobj(n_cell).Color = color_map(round(contour_resolution*contour_mag(n_cell)) == color_index,:); 
        else
            app.gui_plots.contours_gobj(n_cell).Color = app.ops.context_types_all_colors2{tn_all(contour_mag(n_cell))};
        end
    end
end


end