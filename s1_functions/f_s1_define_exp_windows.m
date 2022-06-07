function data = f_s1_define_exp_windows(data, ops)

%% for rest vs grating script define windows
% if ops.processing_type == 1
%     
%     %% will need fixing
%     % Create frame vector indicating the rest phase
%     rest_window = zeros(data.frame_data.num_volumes,1);
%     temp_start = 1;
%     for file_num = 1:ops.num_volt_in_files
%         if frame_data.is_rest(file_num)
%             rest_window(temp_start:(temp_start + data.frame_data.total_frames(file_num) - 1)) = 1;
%         end
%         temp_start = frame_data.total_frames(file_num) + temp_start;
%     end
% 
%     % Define the gratings window
% 
%     fprintf('Defining gratings window...\n');
% 
%     % ask if a gratings window needs to be defined manually
%     figure;
%     plot(data.volt_data_binned_superpos(:,1));
%     hold on;
%     plot(data.ave_trace_superpos);
%     plot(data.volt_data_binned_superpos(:,ops.align_to_channel));
%     title('Do you need to manually define gratings windows? [Y/N]')
%     legend('DAQ voltage trace', 'Mean Ca trace from video', 'LED');
% 
%     define_grating = ask_yes_no_fig();
%     close;
%     if define_grating == 1
%         % initialize
%         grating_window = zeros(data.frame_data.num_volumes,1);
% 
%         for file_num = 1:ops.num_volt_in_files
%             if data.frame_data.is_rest(file_num) == 0
%                 plot_start = sum(frame_data.total_frames(1:file_num)) - frame_data.total_frames(file_num)+1;
%                 plot_end = sum(frame_data.total_frames(1:file_num));
% 
%                 figure;
%                 plot(plot_start:plot_end,volt_data_binned(plot_start:plot_end,1));
%                 hold on;
%                 plot(plot_start:plot_end,norm_redchan1(plot_start:plot_end,1)*3);
%                 plot(plot_start:plot_end,volt_data_binned(plot_start:plot_end,2));
%                 title(sprintf('Select the window for %s (2 clicks)', files_volt_in{file_num}));
%                 legend('DAQ voltage trace', 'Mean Ca trace from video', 'LED');
% 
%                 [x,~]=ginput(2);
% 
%                 if x(1) < plot_start
%                     x(1) = plot_start;
%                 end
%                 if x(2) > plot_end
%                     x(2) = plot_end;
%                 end
%                 x = round(x);
%                 close;
% 
%                 % fill the grating window
%                 grating_window(x(1):x(2)) = 1;            
%             end
%         end
% 
%     else
%         grating_window = 1-rest_window;
%     end
%     close;
%     data.grating_window = grating_window;
%     data.rest_window = rest_window;
% end

%% for auditory script define phases 1 2 3 for control MMN and flipMMN

%if ops.processing_type == 2 || ops.processing_type == 3
exp_phase = zeros(data.frame_data.num_frames_all,1);

if ops.exp_win_selection == 2
    if ~isfield(data, 'exp_window_manual')
        figure;
        hold on;
        plot(data.volt_data_binned_superpos(:,1));
        plot(data.ave_trace_superpos(:,1));
        plot(data.volt_data_binned_superpos(:,2));
        axis tight;
        %plot(data.stim_times_volt{1});
        title('How many exp windows to select?');
        num_win = input('How many exp windows to select?: ');
        title(sprintf('Select %d phases (%d clicks)', num_win, 2*num_win));
        legend('DAQ voltage trace', 'Mean Ca trace from video', 'LED'); % , 'stim times'
        [phase_bounds,~]=ginput(2*num_win);
        close;
        
        phase_bounds = reshape(round(phase_bounds), 2, []);
        phase_bounds = min(phase_bounds, data.frame_data.num_frames_all);

        phase_onset = phase_bounds(1,:);
        phase_offset = phase_bounds(2,:);

        exp_window_manual.phase_onset = phase_onset;
        exp_window_manual.phase_offset = phase_offset;
        if exist(ops.file_save_path_full_processing_params, 'file')
            save(ops.file_save_path_full_processing_params, 'exp_window_manual', '-append')
        else
            save(ops.file_save_path_full_processing_params, 'exp_window_manual')
        end
        data.exp_window_manual = exp_window_manual;
    end

    exp_window = data.exp_window_manual;
elseif ops.exp_win_selection == 1

    %max(data.volt_data_binned_superpos(:,ops.align_to_channel))/std(data.volt_data_binned_superpos(:,ops.align_to_channel));

    align_data = data.volt_data_binned_superpos(:,ops.align_to_channel);
    if sum(align_data == 1) > 3

        pulse_buff = 1; % in sec
        pulse_buff_frames = round(pulse_buff*1000/data.frame_data.frame_period_ave);
        [pulses_onset, pulses_offset] = f_get_pulse_times(data.volt_data_binned_superpos(:,ops.align_to_channel), 0.8);

        phase_onset = pulses_offset + pulse_buff_frames;
        phase_onset(end) = [];
        phase_offset = pulses_onset - pulse_buff_frames;
        phase_offset(1) = [];
    else
        phase_onset = 1;
        phase_offset = numel(align_data);
    end

    exp_window.phase_onset = phase_onset;
    exp_window.phase_offset = phase_offset;
elseif ~ops.exp_win_selection
    exp_window.phase_onset = 1;
    exp_window.phase_offset = size(data.volt_data_binned_superpos,1);
end
data.exp_window = exp_window;
for n_phase = 1:numel(data.exp_window.phase_onset)
    exp_phase(data.exp_window.phase_onset(n_phase):data.exp_window.phase_offset(n_phase)) = n_phase;
end
data.exp_phase = exp_phase;

%     if ops.plot_details
%         figure; hold on;
%         plot(data.volt_data_binned_superpos(:,1))
%         plot(data.exp_phase)
%         tite
%     end
    
%end
data.exp_phase_mpl = f_s1_multiplane_split(exp_phase, ops.num_planes);

end