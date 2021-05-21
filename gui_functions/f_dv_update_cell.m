function f_dv_update_cell(app)

n_pl = app.mplSpinner.Value;

num_cells = app.ddata.num_cells_pl{n_pl};
app.CellSpinner.Value = min([app.CellSpinner.Value, num_cells]);
n_cell = app.CellSpinner.Value;
plot_t = app.ddata.proc_data{1}.frame_data.frame_times_mpl{n_pl}/1000;

%%
app.current_cell_raw = app.cdata.raw(n_cell,:);
app.current_cell_C = app.cdata.C(n_cell,:);
app.current_cell_spikes = app.cdata.S(n_cell,:);

if app.RawButton.Value
    app.gui_plots.plot_raw.XData = plot_t;
    app.gui_plots.plot_raw.YData = app.current_cell_raw;
else
    app.gui_plots.plot_raw.XData = 0;
    app.gui_plots.plot_raw.YData = 0;
end
if app.CButton.Value
    app.gui_plots.plot_C.XData = plot_t;
    app.gui_plots.plot_C.YData = app.current_cell_C;
else
    app.gui_plots.plot_C.XData = 0;
    app.gui_plots.plot_C.YData = 0;
end
if app.SpikesButton.Value
    app.gui_plots.plot_spikes.XData = plot_t;
    app.gui_plots.plot_spikes.YData = app.current_cell_spikes*app.ScaleEditField_sp.Value + app.ShiftEditField_sp.Value;
else
    app.gui_plots.plot_spikes.XData = 0;
    app.gui_plots.plot_spikes.YData = 0;
end

if app.StimtimesButton.Value
    if strcmpi(app.trialtypeDropDown.Value, 'all')
        idx_stim = app.ddata.stim_frame_index{n_pl};
    else
        idx_ctx = strcmpi(app.trialtypeDropDown.Value, app.ops.context_types_labels);
        tt = app.ops.context_types_all(idx_ctx);
        idx_stim = app.ddata.stim_frame_index{n_pl}(app.ddata.trial_types{1} == tt);
    end
    
    stim_trace = zeros(numel(plot_t),1);
    stim_trace(idx_stim) = 1;
    
    app.gui_plots.plot_stim_times.XData = plot_t;
    app.gui_plots.plot_stim_times.YData = stim_trace*app.ScaleEditField_st.Value+app.ShiftEditField_st.Value;
else
    app.gui_plots.plot_stim_times.XData = 0;
    app.gui_plots.plot_stim_times.YData = 0;
end

contours_accepted = app.ddata.OA_data{n_pl}.est.contours(app.ddata.OA_data{n_pl}.proc.comp_accepted);
temp_contours = contours_accepted{n_cell};

if isgraphics(app.gui_plots.plot_current_contour)
    delete(app.gui_plots.plot_current_contour);
end

hold(app.UIAxes, 'on');
app.gui_plots.plot_current_contour = plot(app.UIAxes, temp_contours(:,1), temp_contours(:,2), 'color', [0.75, 0, 0.75], 'LineWidth', 2);
hold(app.UIAxes, 'off');

SNR_accepted = app.ddata.OA_data{n_pl}.proc.SNR2_vals(app.ddata.OA_data{n_pl}.proc.comp_accepted);
app.SNREditField.Value = SNR_accepted(n_cell);

%%

if app.UpdatefigsCheckBox.Value
    if isgraphics(app.gui_plots.freq_resp_fig)
        f_dv_plot_freq_resp(app)
    end
    if isgraphics(app.gui_plots.ctx_resp_fig)
        f_dv_plot_ctx_resp(app)
    end
    if isgraphics(app.gui_plots.select_resp_fig)
        f_dv_plot_select_trial(app)
    end
end

end