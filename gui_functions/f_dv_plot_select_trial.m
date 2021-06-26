function f_dv_plot_select_trial(app)

n_pl = app.mplSpinner.Value;
n_cell = app.CellSpinner.Value;

firing_rate = app.current_cell_spikes;
trial_types = app.ddata.trial_types{1};
stim_times = app.ddata.stim_frame_index{n_pl};
mmn_freq = app.ddata.MMN_freq{1};
trig_window = app.working_ops.trial_num_baseline_resp_frames;
plot_t = app.working_ops.trial_window_t;
num_cont_trials = 400;
num_trials = num_cont_trials/app.ddata.proc_data{1}.stim_params.num_freqs;

trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trig_window);
[trial_data_sort_wctx, trial_types_wctx] = f_s3_add_ctx_trials(trial_data_sort, trial_types, mmn_freq, app.ops);

stats1 = app.ddata.stats{n_pl};

if app.ConverttoZCheckBox.Value
    pop_mean_val = stats1.pop_mean_val(n_cell);
    pop_z_factor = stats1.pop_z_factor(n_cell);
else
    pop_mean_val = 0;
    pop_z_factor = 1;
end

trial_data_sort = (trial_data_sort - pop_mean_val)/pop_z_factor;
trial_data_sort_wctx = (trial_data_sort_wctx - pop_mean_val)/pop_z_factor;

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

pop_mean_trace = mean(trial_data_sort(:,:,1:num_cont_trials),3);
pop_sem_trace = std(trial_data_sort(:,:,1:num_cont_trials), [],3)/sqrt(num_trials-1);
%pop_mean_trace = stats1.pop_mean_trace(n_cell,:);
%pop_sem_trace = stats1.pop_sem_trace(n_cell,:);
%stat_window_t = stats1.stat_window_t;
%stat_plot_intsc = logical(logical(sum(stat_window_t'>=plot_t,2)).*logical(sum(stat_window_t'<=plot_t,2)));
cell_is_resp = stats1.cell_is_resp(n_cell,:);

hold on; axis tight;
plot(plot_t, resp_tr, 'color', [.6 .6 .6])
plot(plot_t, pop_mean_trace, 'color', [0.75, 0, 0.75], 'LineWidth', 2);
plot(plot_t, pop_mean_trace+pop_sem_trace*stats1.z_thresh, '--','color', [0.75, 0, 0.75], 'LineWidth', 1); 
plot(plot_t, mean(resp_tr,2), 'color', [0 0 0], 'LineWidth', 2);
if ~strcmpi(app.trialtypeDropDown.Value, 'all')
    if cell_is_resp(n_tt)
        plot(stats1.peak_t_all(n_cell,n_tt), (stats1.peak_val_all(n_cell,n_tt)-pop_mean_val)/pop_z_factor, '*g')
    end
end
if app.ConverttoZCheckBox.Value
    ylabel('Z scores');
else
    ylabel('Normalized response');
end
title(sprintf('Stim %s; Dset %s; Cell %d', app.trialtypeDropDown.Value, app.ddata.experiment{1}, n_cell), 'Interpreter', 'none')

end