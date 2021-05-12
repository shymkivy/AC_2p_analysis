function f_dv_set_contorus(app)

n_pl = app.mplSpinner.Value;

accepted_cells = app.ddata.OA_data{n_pl}.proc.comp_accepted;
num_cells = sum(accepted_cells);

if strcmp(app.ContoursButtonGroup.SelectedObject.Text,'None')
    app.gui_ops.contour_params.visible_set = 0;
    app.gui_ops.contour_params.contour_mag = zeros(num_cells,1);
    app.gui_ops.contour_params.c_lim = [0 0];
    app.gui_ops.contour_params.c_abs_lim = [0, 0];
elseif strcmp(app.ContoursButtonGroup.SelectedObject.Text,'Tuning')
    if strcmpi(app.trialtypeDropDown.Value, 'all')
        peakvals = max(app.ddata.stats{1}{n_pl}.peak_val_all,[],2);
    else
        ctx_idx = strcmpi(app.ops.context_types_labels, app.trialtypeDropDown.Value);
        peakvals = app.ddata.stats{1}{n_pl}.peak_val_all(:,ctx_idx);
    end
    app.gui_ops.contour_params.visible_set = peakvals>0;
    app.gui_ops.contour_params.contour_mag = peakvals; % factor to resize magnitudes to fir color
    app.gui_ops.contour_params.c_abs_lim = [min(peakvals) max(peakvals)];
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


c_lim = app.gui_ops.contour_params.c_lim;
visible_set = app.gui_ops.contour_params.visible_set;

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

% update colors
if numel(visible_set)>1
    for n_cell = 1:num_cells
        app.gui_plots.contours_gobj(n_cell).Visible = visible_set(n_cell);
        app.gui_plots.contours_gobj(n_cell).Color = color_map(round(contour_resolution*contour_mag(n_cell)) == color_index,:); 
    end
else
    for n_cell = 1:num_cells
        app.gui_plots.contours_gobj(n_cell).Visible = visible_set;
        app.gui_plots.contours_gobj(n_cell).Color = color_map(round(contour_resolution*contour_mag(n_cell)) == color_index,:); 
    end
end


end