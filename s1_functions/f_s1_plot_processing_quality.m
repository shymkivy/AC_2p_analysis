function f_s1_plot_processing_quality(data, ops)

if ops.plot_details
    figure; hold on;
    for n_pl = 1:ops.num_planes
        plot(diff(data.frame_data.frame_times_mpl_linear{n_pl}))
    end
    title('Frame periods across experiments');
    ylabel('Time, ms');
    xlabel('Frame number');
end

% plot the aligned traces with phases and stim times
%data.t_frames_all = 1:data.num_frames_all;
figure; hold on;
if ops.processing_type == 1
    plot(data.frame_data.frame_times_linear,data.volt_data_binned_superpos(:,1));
    plot(data.frame_data.frame_times_linear,data.volt_data_binned_superpos(:,2));
    plot(data.frame_data.frame_times_linear, data.ave_trace_superpos);
    plot(data.frame_data.frame_times_linear,data.rest_window);
    plot(data.frame_data.frame_times_linear,data.grating_window);
    legend('Stim voltage', 'LED', 'Mean Ca trace from video', 'Rest phase', 'Grating phase');
elseif ops.processing_type == 2
    plot(data.frame_data.frame_times_linear,data.volt_data_binned_superpos(:,1));
    plot(data.frame_data.frame_times_linear,data.volt_data_binned_superpos(:,2));
    plot(data.frame_data.frame_times_linear, data.ave_trace_superpos);
    plot(data.frame_data.frame_times_mpl_linear{1}, data.stim_times_trace{1});
    plot(data.frame_data.frame_times_linear, data.mmn_phase);
    legend('Stim voltage', 'LED', 'Mean Ca trace from video', 'Stim start times pl1', 'MMN phase');
elseif ops.processing_type == 3
    plot(data.frame_data.frame_times_linear,data.TDR_trace_binned_superpos);
    plot(data.frame_data.frame_times_linear,data.volt_data_binned_superpos(:,2));
    plot(data.frame_data.frame_times_linear, data.ave_trace_superpos);
    plot(data.frame_data.frame_times_mpl_linear{1}, data.stim_times_trace{1});
    plot(data.frame_data.frame_times_linear, data.mmn_phase);
    legend('TDT voltage', 'LED', 'Mean Ca trace from video', 'Stim start times pl1', 'MMN phase');
end
title('Check everything');

% check alignments for each plane
figure;
for n_pl = 1:ops.num_planes
    temp_volt_data = data.volt_data_binned{n_pl}(:,ops.align_to_channel);
    
    baseline_frames = 5;
    stim_frames = 10;
    t_frames = (-baseline_frames+1):stim_frames;
    
    pulse_times_trace = f_s1_get_stim_onsets(temp_volt_data, 0.5);
    pulse_resp = f_get_stim_trig_resp(data.ave_trace{n_pl}', find(pulse_times_trace), [baseline_frames stim_frames]);
    
    subplot(1,ops.num_planes,n_pl); hold on;
    plot(t_frames, squeeze(pulse_resp));
    xlabel('frame number');
    title(sprintf('Plane %d', n_pl));
    line([1 1], [0 1],'Color','red');
end
sgtitle('Pulse triggered responses');


end