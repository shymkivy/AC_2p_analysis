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
num_frames_mpl = data.frame_data.num_frames_mpl(1);
figure; hold on;
plot(frame_times, volt_data(:,1));
plot(frame_times, volt_data(:,2));
plot(frame_times, data.ave_trace_superpos);
plot(frame_times, data.exp_phase);
legend1 = {'Stim voltage', 'LED', 'Ca trace', 'exp phase'};
if sum(strcmpi(ops.paradigm, {'freq_grating'}))
    plot(frame_times, volt_data(:,4));
    legend1 = [legend1, {'TDT auditory'}];
end

num_chan = numel(data.stim_chan_idx);
for n_ch = 1:num_chan
    if data.stim_chan_idx(n_ch)
        scale = (n_ch)/num_chan + 1;
        stim_times_trace = zeros(num_frames_mpl,1);
        stim_times_trace(data.stim_times_frame{n_ch}) = 1;
        stem(data.frame_data.frame_times_mpl{1}, stim_times_trace*scale, '.', 'markersize', 20);
        legend1 = [legend1, {ops.chan_labels{n_ch}}];
    end
end
axis tight;
legend(legend1);
title(sprintf('Check everything; %s', ops.file_core), 'interpreter', 'none');

% check alignments for each plane

%figure; plot(data.volt_data_all)

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
    axis tight;
end
sgtitle(sprintf('Pulse triggered responses; %s', ops.file_core), 'interpreter', 'none');

end