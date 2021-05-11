function f_dv_update_contour_lims(app)

if ~isempty(app.gui_ops.contour_params.c_lim)
    if strcmp(app.ContoursButtonGroup.SelectedObject.Text,'Tuning')
        app.gui_ops.tuning_lim = app.gui_ops.contour_params.c_lim;
    elseif strcmp(app.ContoursButtonGroup.SelectedObject.Text,'SNR')
        app.gui_ops.SNR_lim = app.gui_ops.contour_params.c_lim;
    end
end

end