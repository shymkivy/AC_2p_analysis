function data = f_s1_process_volt(data, ops)

%% se

if ops.processing_type == 3
    if and(isfield(data.stim_params, 'stim_index'),isfield(data.stim_params, 'sig_dt'))
        if ~isfield(data, 'fg_first_stim_onset')
            figure;
            dcm_obj = datacursormode;
            set(dcm_obj,'UpdateFcn', @NewCallback_YS)
            plot(data.volt_data_all_aligned(:,4));
            hold on;
            plot(data.volt_data_all_aligned(:,ops.align_to_channel));
            title('What is the point of first stimulus onset');
            legend('DAQ voltage trace', 'Alignment channel');

            % ask how many points to use for alignment
            fg_first_stim_onset = input('What is the point of first stimulus onset?:');
            close;
            
            data.fg_first_stim_onset = fg_first_stim_onset;
            if exist(ops.file_save_path_full_processing_params, 'file')
                save(ops.file_save_path_full_processing_params, 'fg_first_stim_onset', '-append')
            else
                save(ops.file_save_path_full_processing_params, 'fg_first_stim_onset')
            end
        end
        stim_index = data.stim_params.stim_index;
        sig_dt = data.stim_params.sig_dt;
        
        stim_times = round(stim_index*sig_dt*1000 + data.fg_first_stim_onset);

        data.stim_times_volt = stim_times;
        data.stim_chan = 4;
    else
        error('No stim_index or sig_dr for freq grating stim_data available');
    end 
end

%% get stim times

% find stim times in volt channels
stim_chan_align = [1]; % 1 5 6
thresh = 0.2;

volt_data = data.volt_data_all_aligned;
volt_time = data.volt_time_all;
[T, num_ch] = size(volt_data);

stim_chan_align = stim_chan_align(stim_chan_align<=num_ch);

stim_times_volt = cell(numel(stim_chan_align),1);
for ch_idx = 1:numel(stim_chan_align)
    n_ch = stim_chan_align(ch_idx);
    stim_trace = volt_data(:,n_ch);
    stim_times_trace = false(T, 1);
    for n_t = 2:T
        if and(stim_trace(n_t) >= thresh, stim_trace(n_t-1) < thresh)
            stim_times_trace(n_t) = 1;
        end
    end
    stim_times_volt{ch_idx} = volt_time(stim_times_trace);
    if numel(stim_times_volt{ch_idx})
        fprintf('Extracted %d stimuli from voltage for chan %d\n', numel(stim_times_volt{ch_idx}),n_ch);
    end
end
if isfield(data, 'stim_chan')
    data.stim_chan = [data.stim_chan, stim_chan_align];
    data.stim_times_volt = [data.stim_times_volt, stim_times_volt];
else
    data.stim_chan = stim_chan_align;
    data.stim_times_volt = stim_times_volt;
end

%% bin voltage data into frame bins
volt_data_binned = cell(1,ops.num_planes);
if isfield(data, 'stim_chan')
    stim_times_frame = cell(numel(data.stim_chan), ops.num_planes);
end

num_vol_lin = data.frame_data.num_volumes_linear;
volt_data = data.volt_data_all_aligned;
num_chan = size(volt_data,2);

for n_pl = 1:ops.num_planes
    volt_data_binned{n_pl} = zeros(num_vol_lin, num_chan);

    temp_frame_times = data.frame_data.frame_times_mpl{n_pl};
    frame_bin_width = mean(data.frame_data.frame_period_ave);
    for n_frame = 1:num_vol_lin
        
        frame_start_index = round(temp_frame_times(n_frame)-frame_bin_width)+1;
        frame_end_index = round(temp_frame_times(n_frame));
        
        % averages visual stim voltage trigger
        for n_ch = 1:num_chan
            if n_ch == 3
                volt_data_binned{n_pl}(n_frame, n_ch) = mean(volt_data(frame_start_index:frame_end_index,n_ch));
            else
                volt_data_binned{n_pl}(n_frame, n_ch) = median(volt_data(frame_start_index:frame_end_index,n_ch));
            end
        end
        
    end
    % process the movement channel
    volt_data_binned{n_pl}(:,3) = abs(gradient(volt_data_binned{n_pl}(:, 3)));
    %data.indexed_volt_data{n_pl} = round(volt_data_binned{n_pl}(:,1)/max(volt_data_binned{n_pl}(:,1))*data.stim_params.num_freqs);
    
    if isfield(data, 'stim_chan')
        for n_ch = 1:numel(data.stim_chan)
            stim_times_ms = data.stim_times_volt{n_ch};
            temp_stim_times_frame = zeros(numel(stim_times_ms), 1);
            for n_st = 1:numel(stim_times_ms)
                temp_st_time = stim_times_ms(n_st);
                temp_stim_times_frame(n_st) = find(temp_frame_times>temp_st_time,1);
            end
            stim_times_frame{n_ch, n_pl} = temp_stim_times_frame;
        end
    end
    
end

% combine over multiplane data
data.volt_data_binned_superpos = f_s1_multiplane_combine(volt_data_binned);

% check if stim start times were binned properly
if ops.processing_type == 3
    if numel(stim_times_frame{end,1}) ~= numel(stim_index)
        error('ERROR!!! Stim times binning didnt work right, line 447');
    end
    
end
if isfield(data, 'stim_chan')
    data.stim_times_frame = stim_times_frame;
end

%% process locomotion here

fprintf('Processing locomotion window...\n');
% process locomotion
if ~isfield(data, 'loco_thresh')
    if ~ops.auto_loco_thresh
        figure;
        plot(data.volt_data_binned_superpos(:,3));
        title('absolute sensor disturbance. pick cutoff for declaring frame a "movment" frame');

        [~,loco_thresh]=ginput(1);
        close;
    else
        loco_thresh = ops.auto_loco_thresh;
    end
    data.loco_thresh = loco_thresh;
    if exist(ops.file_save_path_full_processing_params, 'file')
        save(ops.file_save_path_full_processing_params, 'loco_thresh', '-append')
    else
        save(ops.file_save_path_full_processing_params, 'loco_thresh')
    end
end

for n_pl = 1:ops.num_planes
    temp_data = volt_data_binned{n_pl}(:,3);
    temp_data(temp_data>data.loco_thresh) = 1;
    temp_data(temp_data<=data.loco_thresh) = 0;
    volt_data_binned{n_pl}(:,3) = temp_data;
end

data.volt_data_binned = volt_data_binned;
end
