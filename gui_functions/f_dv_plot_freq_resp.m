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

tn_all = 1:10;
[resp_cells, ~, resp_vals, loc1] = f_dv_get_resp_vals_cells(app, stats1, tn_all, [], 'Resp split');

if app.ConverttoZCheckBox.Value
    st_mean_mean = stats1.stat_trials_mean_mean(n_cell);
    st_mean_sem = stats1.stat_trials_mean_sem(n_cell);
else
    st_mean_mean = 0;
    st_mean_sem = 1;
end

stat_window_idx = and(stats1.stat_window_t>=plot_t(1), stats1.stat_window_t<=plot_t(end));
stat_window_t = stats1.stat_window_t(stat_window_idx);

trial_ave_trace = stats1.stat_trials_mean(n_cell,stat_window_idx);
trial_sem_trace = stats1.stat_trials_sem(n_cell,stat_window_idx);

trial_data_sort = (trial_data_sort - st_mean_mean)/st_mean_sem;
trial_ave_trace = (trial_ave_trace-st_mean_mean)/st_mean_sem;%mean(trial_data_sort(:,:,1:num_cont_trials),3);
trial_sem_trace = trial_sem_trace/st_mean_sem;%std(trial_data_sort(:,:,1:num_cont_trials), [],3)/sqrt(num_trials-1);

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

%trial_ave_trace = stats1.stat_trials_mean(n_cell,:);
%trial_sem_trace = stats1.stat_trials_sem(n_cell,:);
%stat_window_t = stats1.stat_window_t;
%stat_plot_intsc = logical(logical(sum(stat_window_t'>=plot_t,2)).*logical(sum(stat_window_t'<=plot_t,2)));
cell_is_resp = resp_cells(n_cell,:);
for n_tr = 1:10
    subplot(2,5,n_tr); 
    hold on; axis tight; ylim([y_lim_min, y_lim_max]);
    if numel(resp_freq{n_tr})
        if app.IndividualtrialsCheckBox.Value
            plot(plot_t, resp_freq{n_tr}, 'color', [.6 .6 .6]);
        end
        plot(stat_window_t, trial_ave_trace, 'color', [0.75, 0, 0.75], 'LineWidth', 2);
        plot(stat_window_t, trial_ave_trace+trial_sem_trace*stats1.stat_params.z_thresh, '--','color', [0.75, 0, 0.75], 'LineWidth', 1); 
        plot(plot_t, mean(resp_freq{n_tr},2), 'color', [0 0 0], 'LineWidth', 2);
        if cell_is_resp(n_tr)
            if numel(loc1) > 1
                loc2 = loc1(n_cell,n_tr);
            else
                loc2 = loc1;
            end
            plot(loc2, (resp_vals(n_cell,n_tr)-st_mean_mean)/st_mean_sem, '*g')
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
    else
        title(sprintf('Freq %d; no data', n_tr))
    end
end
sgtitle(sprintf('Dset %s; Cell %d', app.ddata.dset_name_full{1}, n_cell), 'Interpreter', 'none');

end