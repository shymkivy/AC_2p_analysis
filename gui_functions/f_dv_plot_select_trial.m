function f_dv_plot_select_trial(app)


n_pl = app.mplSpinner.Value;
n_cell = app.CellSpinner.Value;

firing_rate = app.current_cell_spikes;
trial_types = app.ddata.trial_types{1};
stim_times = app.ddata.stim_frame_index{n_pl};
trig_window = app.working_ops.trial_num_baseline_resp_frames;
plot_t = app.working_ops.trial_window_t;

if strcmpi(app.trialtypeDropDown.Value, 'all')
    idx_stim = stim_times;
else
    idx_ctx = strcmpi(app.trialtypeDropDown.Value, app.ops.context_types_labels);
    tt = app.ops.context_types_all(idx_ctx);
    idx_stim = stim_times(trial_types == tt);
end

temp_resp = f_get_stim_trig_resp(firing_rate, idx_stim, trig_window);
resp_tr = squeeze(temp_resp);

if app.NewplotsCheckBox.Value
    app.gui_plots.select_resp_fig = figure;
else
    if isgraphics(app.gui_plots.select_resp_fig)
        figure(app.gui_plots.select_resp_fig);
        clf(app.gui_plots.select_resp_fig);
    else
        app.gui_plots.select_resp_fig = figure;
    end
end

hold on; axis tight;
plot(plot_t, resp_tr, 'color', [.6 .6 .6])
plot(plot_t, app.stats.resp_all_mean, 'color', [0 0 0], 'LineWidth', 2);
plot(plot_t, app.stats.resp_all_mean+app.stats.z_factor*2, '--','color', [0 0 0], 'LineWidth', 1); 
plot(plot_t, mean(resp_tr,2), 'color', [1 0 1], 'LineWidth', 2);
title(sprintf('Stim %s; Dset %s; Cell %d', app.trialtypeDropDown.Value, app.ddata.experiment{1}, n_cell), 'Interpreter', 'none')


end