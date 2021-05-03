function data = f_s1_process_volt(data, ops)

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
        stim_times_trace_volt = zeros(size(data.volt_data_all_aligned,1),1);
        stim_times_trace_volt(stim_times) = 1;

        data.stim_times_volt = stim_times;
        data.stim_times_trace_volt = stim_times_trace_volt;

    %     figure;
    %     plot(volt_data_all_aligned(:,1));
    %     hold on;
    %     plot(volt_data_all_aligned(:,2));
    %     plot(stim_times_trace);
    %     title('Check if stim onset times are good');
    %     legend('DAQ voltage trace', 'Alignment channel', 'Stim onset times');
    else
        error('No stim_index or sig_dr for freq grating stim_data available');
        return;
    end 
end

%% bin voltage data into frame bins

volt_data_binned = cell(1,ops.num_planes);
if ops.processing_type == 3
   stim_times_trace = cell(1, ops.num_planes);
   TDR_trace_binned = cell(1, ops.num_planes);
end
data.indexed_volt_data = cell(1,ops.num_planes);

for n_pl = 1:ops.num_planes
    volt_data_binned{n_pl} = zeros(data.frame_data.num_volumes_linear, 3);
    if ops.processing_type == 3
        stim_times_trace{n_pl} = zeros(data.frame_data.num_volumes_linear, 1);
        TDR_trace_binned{n_pl} = zeros(data.frame_data.num_volumes_linear, 1);
    end
    
    temp_frame_times = data.frame_data.frame_times_mpl{n_pl};
    frame_bin_width = mean(data.frame_data.frame_period_ave);
    for n_frame=1:data.frame_data.num_volumes_linear
        
        frame_start_index = round(temp_frame_times(n_frame)-frame_bin_width)+1;
        frame_end_index = round(temp_frame_times(n_frame));

        % averages visual stim voltage trigger
        volt_data_binned{n_pl}(n_frame, 1)=median(data.volt_data_all_aligned(frame_start_index:frame_end_index,1),1);
        % averages LED voltage trigger 
        volt_data_binned{n_pl}(n_frame, 2)=median(data.volt_data_all_aligned(frame_start_index:frame_end_index,2),1);
        % takes the absolute val of first derivative
        volt_data_binned{n_pl}(n_frame, 3)=mean(data.volt_data_all_aligned(frame_start_index:frame_end_index,3),1);

        if ops.processing_type == 3
           stim_times_trace{n_pl}(n_frame) = max(stim_times_trace_volt(frame_start_index:frame_end_index));
           TDR_trace_binned{n_pl}(n_frame) = mean(data.volt_data_all_aligned(frame_start_index:frame_end_index,4),1);
        end

    end
    % process the movement channel
    volt_data_binned{n_pl}(:, 3) = abs(gradient(volt_data_binned{n_pl}(:, 3)));
    data.indexed_volt_data{n_pl} = round(volt_data_binned{n_pl}(:,1)/max(volt_data_binned{n_pl}(:,1))*data.stim_params.num_freqs);
end

% combine over multiplane data
data.volt_data_binned_superpos = f_s1_multiplane_combine(volt_data_binned);

% check if stim start times were binned properly
if ops.processing_type == 3
    if sum(stim_times_trace{1}) ~= numel(stim_index)
        error('ERROR!!! Stim times binning didnt work right, line 447');
    end
    data.stim_times_trace = stim_times_trace;
    data.TDR_trace_binned = TDR_trace_binned;
    data.TDR_trace_binned_superpos = f_s1_multiplane_combine(TDR_trace_binned);
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
