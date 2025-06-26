function f_dv_update_contour_lims(app)

contour_val = app.ContoursDropDown.Value;

if ~isempty(app.gui_ops.contour_params.c_lim)
    if strcmpi(contour_val, 'Tuning magnitude')
        app.gui_ops.tuning_lim = app.gui_ops.contour_params.c_lim;
    elseif strcmpi(contour_val, 'SNR')
        app.gui_ops.SNR_lim = app.gui_ops.contour_params.c_lim;
    elseif strcmpi(contour_val, 'skewness')
        app.gui_ops.skew_lim = app.gui_ops.contour_params.c_lim;
    end
end

end