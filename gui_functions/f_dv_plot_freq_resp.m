function f_dv_plot_freq_resp(app)

n_pl = app.mplSpinner.Value;
n_cell = app.CellSpinner.Value;

firing_rate = app.current_cell_spikes;
trial_types = app.ddata.trial_types{1};
stim_times = app.ddata.stim_frame_index{n_pl};

% num_cont_trials = 400;
% num_trials = num_cont_trials/app.ddata.proc_data{1}.stim_params.num_freqs;

trial_window = f_str_to_array(app.plot_BaserespwinEditField.Value);
[plot_t, trial_frames] = f_dv_compute_window_t(trial_window, app.ddata.proc_data{1}.frame_data.volume_period_ave);

trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, trial_frames);

stats1 = app.ddata.stats{n_pl};

if app.ConverttoZCheckBox.Value
    trial_ave_val = stats1.trial_ave_val(n_cell);
    trial_sem_val = stats1.trial_sem_val(n_cell);
else
    trial_ave_val = 0;
    trial_sem_val = 1;
end

stat_window_idx = and(stats1.stat_window_t>=plot_t(1), stats1.stat_window_t<=plot_t(end));
stat_window_t = stats1.stat_window_t(stat_window_idx);

trial_ave_trace = stats1.trial_ave_trace(n_cell,stat_window_idx);
trial_sem_trace = stats1.trial_sem_trace(n_cell,stat_window_idx);

trial_data_sort = (trial_data_sort - trial_ave_val)/trial_sem_val;
trial_ave_trace = (trial_ave_trace-trial_ave_val)/trial_sem_val;%mean(trial_data_sort(:,:,1:num_cont_trials),3);
trial_sem_trace = trial_sem_trace/trial_sem_val;%std(trial_data_sort(:,:,1:num_cont_trials), [],3)/sqrt(num_trials-1);


resp_freq = cell(10,1);
y_lim_max = [0 max(trial_ave_trace+trial_sem_trace*stats1.stat_params.z_thresh)];
y_lim_min = [0 min(trial_ave_trace+trial_sem_trace*stats1.stat_params.z_thresh)];
for n_tr = 1:10
    temp_resp = trial_data_sort(:,:,trial_types == n_tr);
    resp_freq{n_tr} = squeeze(temp_resp);
    if app.IndividualtrialsCheckBox.Value
        y_lim_max = max([y_lim_max max(resp_freq{n_tr}(:))]);
        y_lim_min = min([y_lim_min min(resp_freq{n_tr}(:))]);
    else
        y_lim_max = max([y_lim_max max(mean(resp_freq{n_tr},2))]);
        y_lim_min = min([y_lim_min min(mean(resp_freq{n_tr},2))]);
    end
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


%app.ColoredbytimeCheckBox.Value

%trial_ave_trace = stats1.trial_ave_trace(n_cell,:);
%trial_sem_trace = stats1.trial_sem_trace(n_cell,:);
%stat_window_t = stats1.stat_window_t;
%stat_plot_intsc = logical(logical(sum(stat_window_t'>=plot_t,2)).*logical(sum(stat_window_t'<=plot_t,2)));
cell_is_resp = stats1.cell_is_resp(n_cell,:);
for n_tr = 1:10
    subplot(2,5,n_tr); 
    hold on; axis tight; ylim([y_lim_min, y_lim_max]);
    if app.IndividualtrialsCheckBox.Value
        plot(plot_t, resp_freq{n_tr}, 'color', [.6 .6 .6]);
    end
    plot(stat_window_t, trial_ave_trace, 'color', [0.75, 0, 0.75], 'LineWidth', 2);
    plot(stat_window_t, trial_ave_trace+trial_sem_trace*stats1.stat_params.z_thresh, '--','color', [0.75, 0, 0.75], 'LineWidth', 1); 
    plot(plot_t, mean(resp_freq{n_tr},2), 'color', [0 0 0], 'LineWidth', 2);
    if cell_is_resp(n_tr)
        plot(stats1.peak_t_all(n_cell,n_tr), (stats1.peak_val_all(n_cell,n_tr)-trial_ave_val)/trial_sem_val, '*g')
    end
    if rem(n_tr,5) ~= 1
        set(gca,'ytick',[]);
    else
        if app.ConverttoZCheckBox.Value
            ylabel('Z scores');
        else
            ylabel('Normalized response');
        end
    end
    title(sprintf('Freq %d', n_tr))
end
sgtitle(sprintf('Dset %s; Cell %d', app.ddata.dset_name_full{1}, n_cell), 'Interpreter', 'none');

end