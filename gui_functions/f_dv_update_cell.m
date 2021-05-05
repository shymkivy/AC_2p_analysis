function f_dv_update_cell(app)

n_pl = app.mplSpinner.Value;

num_cells = app.ddata.num_cells_pl{n_pl};
app.CellSpinner.Value = min([app.CellSpinner.Value, num_cells]);
n_cell = app.CellSpinner.Value;


time1 = app.ddata.proc_data{1}.frame_data.frame_times_mpl{n_pl}/1000;


if app.RawButton.Value
    app.plot_raw.XData = time1;
    app.plot_raw.YData = app.ddata.traces_raw{n_pl}(n_cell,:);
else
    app.plot_raw.XData = 0;
    app.plot_raw.YData = 0;
end
if app.CButton.Value
    app.plot_C.XData = time1;
    
    C_all = app.ddata.OA_data{n_pl}.proc.deconv.MCMC.C(app.ddata.OA_data{n_pl}.proc.comp_accepted);
    C1 = C_all{n_cell};
    
    cuts_trace = app.ddata.proc_data{1}.file_cuts_params{n_pl}.vid_cuts_trace;
    C_full = zeros(numel(cuts_trace),1);
    C_full(logical(cuts_trace)) = C1;
    
    app.plot_C.YData = C_full;
else
    app.plot_C.XData = 0;
    app.plot_C.YData = 0;
end
if app.FiringRateButton.Value
    app.plot_fr.XData = time1;
    app.plot_fr.YData = app.ddata.firing_rate{n_pl}(n_cell,:)*app.ScaleEditField_fr.Value+app.ShiftEditField_fr.Value;
else
    app.plot_fr.XData = 0;
    app.plot_fr.YData = 0;
end
if app.SmoothFRButton.Value
    app.plot_fr_sm.XData = time1;
    app.plot_fr_sm.YData = app.ddata.firing_rate_smooth{n_pl}(n_cell,:)*app.ScaleEditField_frsm.Value+app.ShiftEditField_frsm.Value;
else
    app.plot_fr_sm.XData = 0;
    app.plot_fr_sm.YData = 0;
end
if app.StimtimesButton.Value
    if strcmpi(app.trialtypeDropDown.Value, 'all')
        idx_stim = app.ddata.stim_frame_index{n_pl};
    else
        idx_ctx = strcmpi(app.trialtypeDropDown.Value, app.ops.context_types_labels);
        tt = app.ops.context_types_all(idx_ctx);
        idx_stim = app.ddata.stim_frame_index{n_pl}(app.ddata.trial_types{1} == tt);
    end
    
    stim_trace = zeros(numel(time1),1);
    stim_trace(idx_stim) = 1;
    
    app.plot_stim_times.XData = time1;
    app.plot_stim_times.YData = stim_trace*app.ScaleEditField_st.Value+app.ShiftEditField_st.Value;
else
    app.plot_stim_times.XData = 0;
    app.plot_stim_times.YData = 0;
end

contours_accepted = app.ddata.OA_data{n_pl}.est.contours(app.ddata.OA_data{n_pl}.proc.comp_accepted);
temp_contours = contours_accepted{n_cell};

if isgraphics(app.plot_current_contour)
    delete(app.plot_current_contour);
end

hold(app.UIAxes, 'on');
app.plot_current_contour = plot(app.UIAxes, temp_contours(:,1), temp_contours(:,2), 'color', [0.75, 0, 0.75], 'LineWidth', 2);
hold(app.UIAxes, 'off');

SNR_accepted = app.ddata.OA_data{n_pl}.proc.SNR2_vals(app.ddata.OA_data{n_pl}.proc.comp_accepted);
app.SNREditField.Value = SNR_accepted(n_cell);

end