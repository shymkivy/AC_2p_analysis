function f_s1_plot_processing_quality(data, ops)

frame_times = data.frame_data.frame_times_linear;
frame_times_mpl = data.frame_data.frame_times_mpl_linear;
volt_data = data.volt_data_binned_superpos;

if ops.plot_details
    figure; hold on;
    for n_pl = 1:ops.num_planes
        plot(diff(frame_times_mpl{n_pl}))
    end
    title('Frame periods across experiments');
    ylabel('Time, ms');
    xlabel('Frame number');
end

% plot the aligned traces with phases and stim times
%data.t_frames_all = 1:data.num_frames_all;
figure; hold on;
if ops.processing_type == 1
    plot(frame_times, volt_data(:,1));
    plot(frame_times, volt_data(:,2));
    plot(frame_times, data.ave_trace_superpos);
    plot(frame_times, data.rest_window);
    plot(frame_times, data.grating_window);
    legend('Stim voltage', 'LED', 'Mean Ca trace from video', 'Rest phase', 'Grating phase');
else
    
    legend1 = {'LED', 'Mean Ca trace from video', 'MMN phase'};
    if ops.processing_type == 2
        plot(frame_times, volt_data(:,1));
        legend1 = ['Stim voltage', legend1];
    elseif ops.processing_type == 3
        plot(frame_times, volt_data(:,4));
        legend1 = ['TDT voltage', legend1];
    end
    plot(frame_times, volt_data(:,2));
    plot(frame_times, data.ave_trace_superpos);
    plot(frame_times, data.mmn_phase);
    if numel(data.stim_times_frame{1,1})
        stim_times_trace = zeros(numel(frame_times_mpl{1}),1);
        stim_times_trace(data.stim_times_frame{1,1}) = 1;
        plot(frame_times_mpl{1}, stim_times_trace);
        legend1 = [legend1, 'Stim start times pl1'];
    end
    legend(legend1);
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
    pulse_times1 = find(pulse_times_trace);
    if numel(pulse_times1)
        pulse_resp = f_get_stim_trig_resp(data.ave_trace{n_pl}', pulse_times1, [baseline_frames stim_frames]);

        subplot(1,ops.num_planes,n_pl); hold on;
        plot(t_frames, squeeze(pulse_resp));
        xlabel('frame number');
        title(sprintf('Plane %d', n_pl));
        line([1 1], [0 1],'Color','red');
    end
end
sgtitle('Pulse triggered responses');


end