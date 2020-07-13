function data = f_s1_apply_alignment(data, ops)

volt_data_all_segmented_aligned = cell(ops.num_volt_in_files,1);
scaling_factor = data.alignment.scaling_factor;
shift = data.alignment.shift;

for n_file = 1:ops.num_volt_in_files

    % load foltage data
    dat_proc = data.volt_data_all_segmented{n_file};
       
    % align voltage trace here
    shifted_scaled_dat_proc = align_volt_by_scale_shift2(dat_proc, scaling_factor(n_file), shift(n_file));
    dat_proc_len = size(shifted_scaled_dat_proc,1);
    frame_end_time = round(data.frame_data_XML.frame_times_raw{n_file}(end));
    num_chan = size(shifted_scaled_dat_proc,2);
    
    % fill in or cut voltage to make make length the same with XML frames length
    if frame_end_time > dat_proc_len
        last_vals = shifted_scaled_dat_proc(end,:);
        temp_fill = ones(round(frame_end_time - dat_proc_len), num_chan);
        shifted_scaled_dat_proc_pad = [shifted_scaled_dat_proc; temp_fill.*last_vals];
    else
        shifted_scaled_dat_proc_pad = shifted_scaled_dat_proc(1:frame_end_time,:);
    end
    
    % save new aligned traces
    volt_data_all_segmented_aligned{n_file} = shifted_scaled_dat_proc_pad;
end
data.volt_data_all_segmented_aligned = volt_data_all_segmented_aligned;
data.volt_data_all_aligned = cat(1,volt_data_all_segmented_aligned{:});

if ops.plot_details
    figure; 
    plot(volt_data_all_segmented_aligned{1}(:,1));
    hold on;
    for n_pl = 1:ops.num_planes
        plot(data.frame_data.frame_times_mpl{n_pl}, data.ave_trace{n_pl})
    end
    title('check alignment of all planes');
end

%% concatenate frame times files
frame_times_mpl_linear = cell(1,ops.num_planes);
num_frames_mpl_linear = zeros(1,ops.num_planes);
for n_pl = 1:ops.num_planes
    temp_frame_shift = 0;
    frame_times_temp = cell(ops.num_volt_in_files,1);
    for n_file = 1:ops.num_volt_in_files
        % load frame times
        frame_times_temp{n_file} = data.frame_data.frame_times_mpl{n_file,n_pl} + temp_frame_shift;
        % compute frame shift between trials
        temp_frame_shift = frame_times_temp{n_file}(end) + data.frame_data.volume_period(n_file);
    end
    frame_times_mpl_linear{n_pl} = cat(1,frame_times_temp{:});
    num_frames_mpl_linear(n_pl) = numel(frame_times_mpl_linear{n_pl});
end
data.frame_data.frame_times_linear = f_s1_multiplane_combine(frame_times_mpl_linear);
data.frame_data.frame_times_mpl_linear = frame_times_mpl_linear;
data.frame_data.num_frames_mpl_linear = num_frames_mpl_linear;
data.frame_data.num_volumes_linear = min(num_frames_mpl_linear);
end