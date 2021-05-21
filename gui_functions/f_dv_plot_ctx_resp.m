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
    ctx_plot_list = [18 11:17 20;
                     28 21:27 30]';
else
%     ctx_plot_list = [{mmn_freq(2)} {10+app.ops.redundent_pool_trials} {20};...
%                      {mmn_freq(1)} {20+app.ops.redundent_pool_trials} {30}]';
    ctx_plot_list = [18, 19, 20; ...
                     28, 29, 30]';  
end

trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trig_window);
[trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, app.ops);

[n,m] = size(ctx_plot_list);

resp_ctx = cell(m*n,1);
y_lim_max = 0;
y_lim_min = 0;
for n_stim = 1:(m*n)
%     if iscell(ctx_plot_list(n_stim))
%         tt = app.ops.context_types_all(ctx_plot_list{n_stim})';
%         tt_idx = logical(sum(trial_types_wctx == tt,2));
%     else
    tt = app.ops.context_types_all(ctx_plot_list(n_stim));
    tt_idx = trial_types_wctx == tt;

    temp_resp = trial_data_sort_wctx(:,:,tt_idx);
    resp_ctx{n_stim} = squeeze(temp_resp);
    y_lim_max = max([y_lim_max max(resp_ctx{n_stim}(:))]);
    y_lim_min = min([y_lim_min min(resp_ctx{n_stim}(:))]);
end

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
pop_mean = app.ddata.stats{1}{n_pl}.pop_mean{n_cell};
pop_z_factor = app.ddata.stats{1}{n_pl}.pop_z_factor{n_cell};
stat_window_t = app.ddata.stats{1}{n_pl}.stat_window_t;
stat_plot_intsc = logical(logical(sum(stat_window_t'>=plot_t,2)).*logical(sum(stat_window_t'<=plot_t,2)));
cell_is_resp = app.ddata.stats{1}{n_pl}.cell_is_resp(n_cell,:);
for n_stim = 1:(m*n)
    n_tr = ctx_plot_list(n_stim);
    color2 = app.ops.context_types_all_colors2{n_tr};
    title2 = app.ops.context_types_labels{n_tr};
    subplot(m,n,n_stim); hold on; axis tight; ylim([y_lim_min, y_lim_max]);
    plot(plot_t, resp_ctx{n_stim}, 'color', [.6 .6 .6])
    plot(stat_window_t(stat_plot_intsc), pop_mean(stat_plot_intsc), 'color', [0 1 0], 'LineWidth', 2);
    plot(stat_window_t(stat_plot_intsc), pop_mean(stat_plot_intsc)+pop_z_factor(stat_plot_intsc)*app.ddata.stats{1}{n_pl}.z_thresh, '--','color', [0 1 0], 'LineWidth', 1); 
    plot(stat_window_t(stat_plot_intsc), mean(resp_ctx{n_stim},2), 'color', color2, 'LineWidth', 2);
    if cell_is_resp(n_tr)
        plot(app.ddata.stats{1}{n_pl}.peak_t_all(n_cell,n_tr), app.ddata.stats{1}{n_pl}.peak_val_all(n_cell,n_tr), '*g')
    end
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