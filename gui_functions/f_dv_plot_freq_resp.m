function f_dv_plot_freq_resp(app)

n_pl = app.mplSpinner.Value;
n_cell = app.CellSpinner.Value;

firing_rate = app.current_cell_spikes;
trial_types = app.ddata.trial_types{1};
stim_times = app.ddata.stim_frame_index{n_pl};
trig_window = app.working_ops.trial_num_baseline_resp_frames;
plot_t = app.working_ops.trial_window_t;

resp_freq = cell(10,1);
y_lim_max = 0;
y_lim_min = 0;
for n_tr = 1:10
    temp_resp = f_get_stim_trig_resp(firing_rate, stim_times(trial_types == n_tr), trig_window);
    resp_freq{n_tr} = squeeze(temp_resp);
    y_lim_max = max([y_lim_max max(resp_freq{n_tr}(:))]);
    y_lim_min = min([y_lim_min min(resp_freq{n_tr}(:))]);
end
if y_lim_max == 0
    y_lim_max = 1;
end

if app.NewplotsCheckBox.Value
    app.gui_plots.freq_resp_fig = figure;
else
    if isgraphics(app.gui_plots.freq_resp_fig)
        figure(app.gui_plots.freq_resp_fig);
        clf(app.gui_plots.freq_resp_fig);
    else
        app.gui_plots.freq_resp_fig = figure;
    end
end
pop_mean = app.ddata.stats.pop_mean{n_cell};
pop_z_factor = app.ddata.stats.pop_z_factor{n_cell};
stat_window_t = app.ddata.stats.stat_window_t;
stat_plot_intsc = logical(logical(sum(stat_window_t'>plot_t,2)).*logical(sum(stat_window_t'<plot_t,2)));
cell_is_resp = app.ddata.stats.cell_is_resp(n_cell,:);
for n_tr = 1:10
    subplot(2,5,n_tr); 
    hold on; axis tight; ylim([y_lim_min, y_lim_max]);
    plot(plot_t, resp_freq{n_tr}, 'color', [.6 .6 .6])
    plot(stat_window_t(stat_plot_intsc), pop_mean(stat_plot_intsc), 'color', [0 0 0], 'LineWidth', 2);
    plot(stat_window_t(stat_plot_intsc), pop_mean(stat_plot_intsc)+pop_z_factor(stat_plot_intsc)*2, '--','color', [0 0 0], 'LineWidth', 1); 
    plot(plot_t, mean(resp_freq{n_tr},2), 'color', [1 0 1], 'LineWidth', 2);
    if cell_is_resp(n_tr)
        plot(app.ddata.stats.peak_t_all(n_cell,n_tr), app.ddata.stats.peak_val_all(n_cell,n_tr), '*g')
    end
    if rem(n_tr,5) ~= 1
        set(gca,'ytick',[])
    end
    title(sprintf('Freq %d', n_tr))
end
sgtitle(sprintf('Dset %s; Cell %d', app.ddata.experiment{1}, n_cell), 'Interpreter', 'none')
end