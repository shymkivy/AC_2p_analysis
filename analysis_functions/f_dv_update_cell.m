function f_dv_update_cell(app)

params = f_dv_gather_params(app);
params.ddata = app.ddata;

n_pl = params.n_pl;

cdata = f_dv_compute_cdata(app.ddata, params);
app.cdata{n_pl} = cdata;

num_cells = cdata.num_cells;
app.CellSpinner.Value = min([app.CellSpinner.Value, num_cells]);
n_cell = app.CellSpinner.Value;
plot_t = app.ddata.proc_data{1}.frame_data.frame_times_mpl{n_pl}/1000;
accepted_cells = cdata.accepted_cells;

%%
app.current_cell_raw = cdata.raw(n_cell,:);
app.current_cell_C = cdata.C(n_cell,:);
app.current_cell_spikes = cdata.S_sm(n_cell,:);

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
    tn_all = f_dv_get_trial_number(params);
    tt = app.ops.context_types_all(tn_all)';
    if isempty(app.ddata.stim_frame_index{n_pl}) || isempty(app.ddata.trial_types{1})
        idx_stim = [];
    else
        idx_stim = app.ddata.stim_frame_index{n_pl}(logical(sum(app.ddata.trial_types{1} == tt,2)));
    end
    
    stim_trace = zeros(numel(plot_t),1);
    stim_trace(idx_stim) = 1;
    
    app.gui_plots.plot_stim_times.XData = plot_t;
    app.gui_plots.plot_stim_times.YData = stim_trace*app.ScaleEditField_st.Value+app.ShiftEditField_st.Value;
else
    app.gui_plots.plot_stim_times.XData = 0;
    app.gui_plots.plot_stim_times.YData = 0;
end

contours_accepted = app.ddata.OA_data{n_pl}.est.contours(accepted_cells);
temp_contours = contours_accepted{n_cell};

if isgraphics(app.gui_plots.plot_current_contour)
    delete(app.gui_plots.plot_current_contour);
end

hold(app.UIAxes, 'on');
app.gui_plots.plot_current_contour = plot(app.UIAxes, temp_contours(:,2), temp_contours(:,1), 'color', [0.75, 0, 0.75], 'LineWidth', 2);
hold(app.UIAxes, 'off');

SNR_accepted = app.ddata.OA_data{n_pl}.proc.SNR2_vals(accepted_cells);
app.SNREditField.Value = SNR_accepted(n_cell);

%%
if app.UpdatefigsCheckBox.Value
    if isgraphics(app.gui_plots.freq_resp_fig)
        f_dv_plot_freq_resp(app);
    end
    if isgraphics(app.gui_plots.ctx_resp_fig)
        f_dv_plot_ctx_resp(app);
    end
    if isgraphics(app.gui_plots.select_resp_fig)
        f_dv_plot_select_trial(app);
    end
    if isgraphics(app.gui_plots.tuning_fig)
        f_dv_plot_tuning(app);
    end
end

end