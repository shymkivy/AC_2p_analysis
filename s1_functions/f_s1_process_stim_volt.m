function data = f_s1_process_stim_volt(data, ops)

% if ~isfield(data, 'stim_times_trace')
%     
%     stim_times_trace = cell(1,ops.num_planes);
%     for n_pl = 1:ops.num_planes
%         % Extract stim onset time
%         stim_cutoff = 1/2;
%         stim_times_trace{n_pl} = f_s1_get_stim_onsets(data.indexed_volt_data{n_pl}.*logical(data.mmn_phase_mpl{n_pl}),stim_cutoff);
%         
%         fprintf('Extracted %d stimuli from voltage file for plane %d\n', sum(stim_times_trace{n_pl}), n_pl);
% %         figure; hold on;
% %         plot(indexed_volt_data{n_pl});
% %         plot(stim_times_trace{n_pl});
%     end
%     data.stim_times_trace = stim_times_trace;
%     
% end

%%
% data.stim_frame_index = cell(1,ops.num_planes);
% data.num_trials_extracted = zeros(1,ops.num_planes);
% for n_pl = 1:ops.num_planes
%     data.stim_frame_index{n_pl} = find(data.stim_times_trace{n_pl});
%     data.num_trials_extracted(n_pl) = sum(data.stim_times_trace{n_pl});
% end

%%
if ~isfield(data, 'cont_trials_seq')
    % normalize the trigger voltage to index
    stim_ch_idx = (data.stim_chan == 1);
    stim_frames = data.stim_times_frame{stim_ch_idx,1};
    if numel(stim_frames)
        mmn_phase = data.mmn_phase_mpl{1};
        volt_data = data.volt_data_binned{1}(:,1);

        num_freqs = data.stim_params.num_freqs;
        max_cont_volt = max((mmn_phase == 1) .* (volt_data));

        stim_phase_frames = mmn_phase(stim_frames);

        num_cont = sum(stim_phase_frames == 1);
        num_mmn = sum(stim_phase_frames == 2);

        cont_trials_seq = zeros(num_cont, 1);
        MMN_trials_seq = zeros(num_mmn, 2);

        stim_dur_frames = floor(data.stim_params.stim_duration/data.frame_data.volume_period_ave*1000);

        % Extract stimulus direction
        %trial_type_directions = zeros(numel(stim_frame_index),1);
        for n_trial = 1:numel(stim_frames)
            % do a median over frames to prevent some unlikely noise
            temp_frame = stim_frames(n_trial);
            if stim_phase_frames(n_trial) == 1
                cont_trials_seq(n_trial) = round(median(volt_data(temp_frame:(temp_frame+stim_dur_frames-1)))/max_cont_volt*num_freqs);
            elseif stim_phase_frames(n_trial) > 1
                MMN_trials_seq(n_trial - num_cont) = round(median(volt_data(temp_frame:temp_frame+stim_dur_frames)));
            end
        end
    end
else
    cont_trials_seq = data.cont_trials_seq;
    MMN_trials_seq = data.MMN_trials_seq;
end

%% process redundents and deviants
% MMN redundants are 101 - 140, deviant is 170. 
% flipMMN redundantts are 201 - 240, deviant is 270
% which stim patter did we use. row 1 is MMN, row 2 is flip
if exist('cont_trials_seq', 'var')
    if sum(max(MMN_trials_seq)> 2)
        MMN_trials_seq = round(2*MMN_trials_seq/(max(MMN_trials_seq(:))));
    end

    MMN_trials_seq_indexed = zeros(size(MMN_trials_seq));

    stim_params.red_dev_stim_patten = [mode(MMN_trials_seq(:,1)), 3-mode(MMN_trials_seq(:,1));
                                       mode(MMN_trials_seq(:,2)), 3-mode(MMN_trials_seq(:,2))];

    MMN_trials_seq_indexed(MMN_trials_seq(:,1) == stim_params.red_dev_stim_patten(1,2),1) = 170;
    MMN_trials_seq_indexed(MMN_trials_seq(:,2) == stim_params.red_dev_stim_patten(2,2),2) = 270;

    % mark MMN and flipMMN
    for n_flip = 1:2
        temp_seq = 0;
        for n_trial = 1:size(MMN_trials_seq,1)
            if  MMN_trials_seq(n_trial,n_flip) == stim_params.red_dev_stim_patten(n_flip,1)
                MMN_trials_seq_indexed(n_trial,n_flip) = 1 + n_flip*100 + temp_seq;
                temp_seq = temp_seq + 1;
            else
                temp_seq = 0;
            end
        end
    end

    data.cont_trials_seq = cont_trials_seq;
    data.MMN_trials_seq = MMN_trials_seq_indexed;
    data.trial_types = [cont_trials_seq; MMN_trials_seq_indexed(:,1); MMN_trials_seq_indexed(:,2)];
end
end