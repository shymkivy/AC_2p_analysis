function f_plot_spike_inference_details(data, plot_cell_num, ops) 
% this plots the trial traces + trial averages for all spike inference
% methods, including raw, smoothed dt/dt, spike inerence, and smoothed
% spike inference


trial_type = data.trial_types;
raw_cell_trace = data.cell_trace(plot_cell_num,:);
smooth_dfdt = data.ddt_cell_dataN(plot_cell_num,:);
spike_inf = data.spike_inf(plot_cell_num,:);
smooth_spike_inf = data.smooth_spike_inf(plot_cell_num,:);

stim_times = data.stim_frame_index;
num_frames_full = size(data.cell_trace,2);
time_window = data.time_stim_window;
trial_ave_win = [sum(time_window<=0) sum(time_window>0)];

temp_data_cell = cell(4,1);
temp_data_cell{1} = f_get_stim_trig_resp(raw_cell_trace, stim_times, trial_ave_win);
temp_data_cell{1} = temp_data_cell{1} - mean(temp_data_cell{1}, 2);
temp_data_cell{2} = f_get_stim_trig_resp(smooth_dfdt, stim_times, trial_ave_win);
temp_data_cell{3} = f_get_stim_trig_resp(spike_inf, stim_times, trial_ave_win);
temp_data_cell{4} = f_get_stim_trig_resp(smooth_spike_inf, stim_times, trial_ave_win);


plot_titles = {'Raw', 'Smoothed df/dt', 'Spike inf', ['Smoothed spike inf; sigma=' num2str(ops.spike_inf_gaus_sigma) 'frame']};
for n_cntxt = 1:numel(ops.context_types_all)
    ctx_trials_indx = (trial_type == ops.context_types_all(n_cntxt));
    figure;
    for ii = 1:4
        [temp_trial_ave, baseline] = f_trial_average(temp_data_cell{ii},trial_type,ops);
        %z_factor = std(reshape(temp_trial_ave, [], 1));     
        z_factor = mean(sqrt((temp_trial_ave(:)).^2));
        temp_data = squeeze(temp_data_cell{ii}(1,ops.win.analysis_window,:))-baseline;
        y_lim_max = max(temp_data(:)/z_factor);
        y_lim_min = min(temp_data(:)/z_factor);
        
        subplot(2,4,ii);
        plot(time_window(ops.win.analysis_window), temp_data(:,ctx_trials_indx)/z_factor, 'Color', [.8 .8 .8]);
        hold on;
        plot(time_window(ops.win.analysis_window), squeeze(temp_trial_ave(1,ops.win.analysis_window,n_cntxt))/z_factor,'k','LineWidth',2);
        plot(time_window(ops.win.analysis_window),ones(1,sum(ops.win.analysis_window))*ops.z_threshold, 'r');
        axis tight;
        ylim([y_lim_min-0.1 y_lim_max]);
        title(plot_titles{ii});
        if ii == 1
            xlabel('Time (sec)');
            ylabel('Z-score');
        end
    end
    subplot(2,4,5:8);
    plot(raw_cell_trace); hold on;
    plot(smooth_dfdt);
    plot(spike_inf);
    plot(smooth_spike_inf);
    stim_times_trace = zeros(1,num_frames_full);
    stim_times_trace(stim_times(ctx_trials_indx)) = 5;
    plot(stim_times_trace);
    legend('Raw trace', 'Smoothed df/dt', 'Spike inf', 'Contex stim times');
    title(['Ctx ' num2str(n_cntxt) '; SNR=' num2str(data.SNR_vals(plot_cell_num))]);
    xlabel('Time (frames)');
    
end
end