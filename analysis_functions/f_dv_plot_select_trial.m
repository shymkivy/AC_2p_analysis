function f_dv_plot_select_trial(app)

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

trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trial_frames);
[trial_data_sort_wctx, trial_types_wctx] = f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, app.ops);

stats1 = app.ddata.stats{n_pl};

if app.ConverttoZCheckBox.Value
    st_mean_mean = stats1.stat_trials_mean_mean(n_cell);
    st_mean_sem = stats1.stat_trials_mean_sem(n_cell);
else
    st_mean_mean = 0;
    st_mean_sem = 1;
end

trial_data_sort = (trial_data_sort - st_mean_mean)/st_mean_sem;
trial_data_sort_wctx = (trial_data_sort_wctx - st_mean_mean)/st_mean_sem;

tn_all = f_dv_get_trial_number(app);
tt_all = app.ops.context_types_all(tn_all);
idx1 = logical(sum(trial_types_wctx == tt_all',2));
temp_resp = trial_data_sort_wctx(:,:,idx1);
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

trial_ave_trace = mean(trial_data_sort(:,:,1:num_cont_trials),3);
trial_sem_trace = std(trial_data_sort(:,:,1:num_cont_trials), [],3)/sqrt(num_trials-1);
%trial_ave_trace = stats1.stat_trials_mean(n_cell,:);
%trial_sem_trace = stats1.stat_trials_sem(n_cell,:);
%stat_window_t = stats1.stat_window_t;
%stat_plot_intsc = logical(logical(sum(stat_window_t'>=plot_t,2)).*logical(sum(stat_window_t'<=plot_t,2)));
cell_is_resp = stats1.peak_resp_cells(n_cell,:);

hold on; axis tight;
plot(plot_t, resp_tr, 'color', [.6 .6 .6])
plot(plot_t, trial_ave_trace, 'color', [0.75, 0, 0.75], 'LineWidth', 2);
plot(plot_t, trial_ave_trace+trial_sem_trace*stats1.stat_params.z_thresh, '--','color', [0.75, 0, 0.75], 'LineWidth', 1); 
plot(plot_t, mean(resp_tr,2), 'color', [0 0 0], 'LineWidth', 2);
if ~sum(strcmpi(app.trialtypeDropDown.Value, {'all', 'Freqs', 'Context'}))
    if cell_is_resp(tn_all)
        plot(stats1.peak_loc(n_cell,tn_all), (stats1.peak_vals(n_cell,tn_all)-st_mean_mean)/st_mean_sem, '*g')
    end
end
if app.ConverttoZCheckBox.Value
    ylabel('Z scores');
else
    ylabel('Normalized response');
end
title(sprintf('Stim %s; Dset %s; Cell %d', app.trialtypeDropDown.Value, app.ddata.dset_name_full{1}, n_cell), 'Interpreter', 'none')

end