function f_dv_plot_ctx_resp(app)


n_pl = app.mplSpinner.Value;
n_cell = app.CellSpinner.Value;

firing_rate = app.current_cell_spikes;
trial_types = app.ddata.trial_types{1};
stim_times = app.ddata.stim_frame_index{n_pl};
mmn_freq = app.ddata.MMN_freq{1};
trig_window = app.working_ops.trial_num_baseline_resp_frames;
plot_t = app.working_ops.trial_window_t;

if app.FullredCheckBox.Value
    ctx_plot_list = [mmn_freq(2) 11:17 20; mmn_freq(1) 21:27 30]';
else
    ctx_plot_list = [{mmn_freq(2)} {10+app.ops.redundent_pool_trials} {20};...
                     {mmn_freq(1)} {20+app.ops.redundent_pool_trials} {30}]';
end

[n,m] = size(ctx_plot_list);

resp_ctx = cell(m*n,1);
y_lim_max = 0;
y_lim_min = 0;
for n_stim = 1:(m*n)
    if iscell(ctx_plot_list(n_stim))
        tt = app.ops.context_types_all(ctx_plot_list{n_stim})';
        tt_idx = logical(sum(trial_types == tt,2));
    else
        tt = app.ops.context_types_all(ctx_plot_list(n_stim));
        tt_idx = trial_types == tt;
    end
    temp_resp = f_get_stim_trig_resp(firing_rate, stim_times(tt_idx), trig_window);
    resp_ctx{n_stim} = squeeze(temp_resp);
    y_lim_max = max([y_lim_max max(resp_ctx{n_stim}(:))]);
    y_lim_min = min([y_lim_min min(resp_ctx{n_stim}(:))]);
end

colors1 = [app.ops.context_colors, app.ops.context_colors];
titles1 = {'Cont 2', 'RedF', 'Dev', 'Cont 1', 'Red', 'DevF'};

if app.NewplotsCheckBox.Value
    app.gui_plots.ctx_resp_fig = figure;
else
    if isgraphics(app.gui_plots.ctx_resp_fig)
        figure(app.gui_plots.ctx_resp_fig);
        clf(app.gui_plots.ctx_resp_fig)
    else
        app.gui_plots.ctx_resp_fig = figure;
    end
end
for n_stim = 1:(m*n)
    if iscell(ctx_plot_list(n_stim))
        color2 = colors1{n_stim};
        title2 = titles1{n_stim};
    else
        color2 = app.ops.context_types_all_colors2{ctx_plot_list(n_stim)};
        title2 = app.ops.context_types_labels{ctx_plot_list(n_stim)};
    end
    subplot(m,n,n_stim); hold on; axis tight; ylim([y_lim_min, y_lim_max]);
    plot(plot_t, resp_ctx{n_stim}, 'color', [.6 .6 .6])
    plot(plot_t, app.stats.resp_all_mean, 'color', [0 1 0], 'LineWidth', 2);
    plot(plot_t, app.stats.resp_all_mean+app.stats.z_factor*2, '--','color', [0 1 0], 'LineWidth', 1); 
    plot(plot_t, mean(resp_ctx{n_stim},2), 'color', color2, 'LineWidth', 2);
    if rem(n_stim,n) ~= 1
        set(gca,'ytick',[])
    end
    if rem(n_stim,n) == 1
        if n_stim > n
            ylabel(app.ops.context_types_labels{mmn_freq(1)})
        else
            ylabel(app.ops.context_types_labels{mmn_freq(2)})
        end
    end
    title(sprintf('%s', title2))
end
sgtitle(sprintf('Dset %s; Cell %d', app.ddata.experiment{1}, n_cell), 'Interpreter', 'none')
end