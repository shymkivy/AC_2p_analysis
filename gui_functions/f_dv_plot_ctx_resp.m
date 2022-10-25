function f_dv_plot_ctx_resp(app)

n_pl = app.mplSpinner.Value;
n_cell = app.CellSpinner.Value;

firing_rate = app.current_cell_spikes;
trial_types = app.ddata.trial_types{1};
stim_times = app.ddata.stim_frame_index{n_pl};
mmn_freq = app.ddata.MMN_freq{1};

trial_window = f_str_to_array(app.plot_BaserespwinEditField.Value);
[plot_t, trial_frames] = f_dv_compute_window_t(trial_window, app.ddata.proc_data{1}.frame_data.volume_period_ave);

num_cont_trials = 400;
num_trials = num_cont_trials/app.ddata.proc_data{1}.stim_params.num_freqs;

if app.FullredCheckBox.Value
    ctx_plot_list = [18 11:17 20;
                     28 21:27 30]';
else
%     ctx_plot_list = [{mmn_freq(2)} {10+app.ops.redundent_pool_trials} {20};...
%                      {mmn_freq(1)} {20+app.ops.redundent_pool_trials} {30}]';
    ctx_plot_list = [18, 19, 20; ...
                     28, 29, 30]';  
end

trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trial_frames);
[trial_data_sort_wctx, trial_types_wctx] =  f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, app.ops);

stats1 = app.ddata.stats{n_pl};

if strcmpi(app.TuningfreatureDropDown.Value, 'peaks')
    resp_cells = stats1.peak_resp_cells;
elseif strcmpi(app.TuningfreatureDropDown.Value, 'onset')
    resp_cells = stats1.onset_resp_cells;
elseif strcmpi(app.TuningfreatureDropDown.Value, 'offset')
    resp_cells = stats1.offset_resp_cells;
end

if app.ConverttoZCheckBox.Value
    st_mean_mean = stats1.stat_trials_mean_mean(n_cell);
    st_mean_sem = stats1.stat_trials_mean_sem(n_cell);
else
    st_mean_mean = 0;
    st_mean_sem = 1;
end

trial_data_sort = (trial_data_sort - st_mean_mean)/st_mean_sem;
trial_data_sort_wctx = (trial_data_sort_wctx - st_mean_mean)/st_mean_sem;

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

trial_ave_trace = mean(trial_data_sort(:,:,1:num_cont_trials),3);
trial_sem_trace = std(trial_data_sort(:,:,1:num_cont_trials), [],3)/sqrt(num_trials-1);
%trial_ave_trace = stats1.stat_trials_mean(n_cell,:);
%trial_sem_trace = stats1.stat_trials_sem(n_cell,:);
%stat_window_t = stats1.stat_window_t;
%stat_plot_intsc = logical(logical(sum(stat_window_t'>=plot_t,2)).*logical(sum(stat_window_t'<=plot_t,2)));
cell_is_resp = resp_cells(n_cell,:);
for n_stim = 1:(m*n)
    n_tr = ctx_plot_list(n_stim);
    color2 = app.ops.context_types_all_colors2{n_tr};
    title2 = app.ops.context_types_labels{n_tr};
    subplot(m,n,n_stim); hold on; axis tight; ylim([y_lim_min, y_lim_max]);
    plot(plot_t, resp_ctx{n_stim}, 'color', [.6 .6 .6])
    plot(plot_t, trial_ave_trace, 'color', [0.75, 0, 0.75], 'LineWidth', 2);
    plot(plot_t, trial_ave_trace+trial_sem_trace*stats1.stat_params.z_thresh, '--','color', [0.75, 0, 0.75], 'LineWidth', 1); 
    plot(plot_t, mean(resp_ctx{n_stim},2), 'color', color2, 'LineWidth', 2);
    if cell_is_resp(n_tr)
        if numel(stats2.loc) > 1
            loc1 = stats2.loc(n_cell,n_tr);
        else
            loc1 = stats2.loc;
        end
        plot(loc1, (stats2.vals(n_cell,n_tr)-st_mean_mean)/st_mean_sem, '*g')
    end
    if rem(n_stim,n) ~= 1
        set(gca,'ytick',[])
    end
    if rem(n_stim,n) == 1
        if n_stim > n
            if app.ConverttoZCheckBox.Value
                ylabel(sprintf('%s; %s',app.ops.context_types_labels{mmn_freq(1)}, 'Z scores'));
            else
                ylabel(sprintf('%s; %s',app.ops.context_types_labels{mmn_freq(1)}, 'norm resp'));
            end
        else
            if app.ConverttoZCheckBox.Value
                ylabel(sprintf('%s; %s',app.ops.context_types_labels{mmn_freq(2)}, 'Z scores'));
            else
                ylabel(sprintf('%s; %s',app.ops.context_types_labels{mmn_freq(2)}, 'norm resp'));
            end
        end
    end
    title(sprintf('%s', title2))
end
sgtitle(sprintf('Dset %s; Cell %d', app.ddata.dset_name_full{1}, n_cell), 'Interpreter', 'none')

end