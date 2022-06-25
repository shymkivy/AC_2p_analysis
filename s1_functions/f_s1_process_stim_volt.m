function data = f_s1_process_stim_volt(data, ops)
%%

% if ~isfield(data, 'stim_times_trace')
%     
%     stim_times_trace = cell(1,ops.num_planes);
%     for n_pl = 1:ops.num_planes
%         % Extract stim onset time
%         stim_cutoff = 1/2;
%         stim_times_trace{n_pl} = f_s1_get_stim_onsets(data.indexed_volt_data{n_pl}.*logical(data.exp_phase_mpl{n_pl}),stim_cutoff);
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

if strcmpi(ops.paradigm, {'freq_grating'})
    trial_types_all = cell(3,1);
    trial_types_all{1} = data.cont_trials_seq;
    trial_types_all{2} = data.MMN_trials_seq(:,1);
    trial_types_all{3} = data.MMN_trials_seq(:,2);
elseif sum(strcmpi(ops.paradigm, {'ammn', 'vmmn', 'mmn', 'ammn_stim', 'behavior', 'cont'}))
    % normalize the trigger voltage to index
    stim_ch_idx = strcmpi(ops.chan_labels, 'stim type');
    stim_frames = data.stim_times_frame{stim_ch_idx,1};
    
    num_phases = numel(data.exp_window.phase_onset);
    
    trial_types_all = cell(num_phases,1);
    num_trials_all = zeros(num_phases,1);
    
    exp_phase = data.exp_phase_mpl{1};
    volt_data = data.volt_data_binned{1}(:,stim_ch_idx);
    num_freqs = data.stim_params.num_freqs;
    num_frames = numel(volt_data);
    
    stim_phase_frames = exp_phase(stim_frames);
    
    stim_dur_frames = floor(data.stim_params.stim_duration/data.frame_data.volume_period_ave*1000);
    
    for n_ph = 1:num_phases
        idx1 = stim_phase_frames == n_ph;
        stim_frames2 = stim_frames(idx1);
        num_trials_all(n_ph) = sum(idx1);
        
        trial_types = zeros(num_trials_all(n_ph),1);
        
        for n_trial = 1:num_trials_all(n_ph)
            temp_frame = stim_frames2(n_trial);
            temp_frame_end = min(temp_frame+stim_dur_frames-1, num_frames);
            trial_types(n_trial) = median(volt_data(temp_frame:temp_frame_end));
        end
        if n_ph == 1
            trial_types_all{n_ph} = round(trial_types/4*num_freqs);
        else
            trial_types_all{n_ph} = round(trial_types);
        end
        
    end
    
    
%     if numel(stim_frames)
%         exp_phase = data.exp_phase_mpl{1};
%         volt_data = data.volt_data_binned{1}(:,1);
% 
%         num_freqs = data.stim_params.num_freqs;
%         max_cont_volt = max((exp_phase == 1) .* (volt_data));
% 
%         stim_phase_frames = exp_phase(stim_frames);
% 
%         num_cont = sum(stim_phase_frames == 1);
%         num_mmn = sum(stim_phase_frames == 2);
% 
%         cont_trials_seq = zeros(num_cont, 1);
%         MMN_trials_seq = zeros(num_mmn, 2);
% 
%         stim_dur_frames = floor(data.stim_params.stim_duration/data.frame_data.volume_period_ave*1000);
% 
%         % Extract stimulus direction
%         %trial_type_directions = zeros(numel(stim_frame_index),1);
%         for n_trial = 1:numel(stim_frames)
%             % do a median over frames to prevent some unlikely noise
%             temp_frame = stim_frames(n_trial);
%             if stim_phase_frames(n_trial) == 1
%                 cont_trials_seq(n_trial) = round(median(volt_data(temp_frame:(temp_frame+stim_dur_frames-1)))/max_cont_volt*num_freqs);
%             elseif stim_phase_frames(n_trial) > 1
%                 MMN_trials_seq(n_trial - num_cont) = round(median(volt_data(temp_frame:temp_frame+stim_dur_frames)));
%             end
%         end
%     end
end

%% process redundents and deviants
% MMN redundants are 101 - 140, deviant is 170. 
% flipMMN redundantts are 201 - 240, deviant is 270
% which stim patter did we use. row 1 is MMN, row 2 is flip

if sum(strcmpi(ops.paradigm, {'ammn', 'vmmn', 'mmn', 'ammn_stim'}))
    
    stim_params.red_dev_stim_patten = zeros(2, 2);
    for n_ph = 2:3
        MMN_trials_seq = trial_types_all{n_ph};
        
        if sum(max(MMN_trials_seq)> 2)
            MMN_trials_seq = round(2*MMN_trials_seq/(max(MMN_trials_seq(:))));
        end
        
        MMN_trials_seq_indexed = zeros(size(MMN_trials_seq));
        
        red_dev_stim_patten = [mode(MMN_trials_seq), 3-mode(MMN_trials_seq)];
        stim_params.red_dev_stim_patten(n_ph-1, :) = red_dev_stim_patten;
                                       
        MMN_trials_seq_indexed(MMN_trials_seq == red_dev_stim_patten(2)) = (n_ph-1)*100+70;
               
        temp_seq = 0;
        for n_trial = 1:numel(MMN_trials_seq)
            if  MMN_trials_seq(n_trial) == stim_params.red_dev_stim_patten(1)
                MMN_trials_seq_indexed(n_trial) = 1 + (n_ph-1)*100 + temp_seq;
                temp_seq = temp_seq + 1;
            else
                temp_seq = 0;
            end
        end
        trial_types_all{n_ph} = MMN_trials_seq_indexed;
    end
end

if exist('trial_types_all', 'var')
    data.trial_types_all = trial_types_all;
    data.trial_types = cat(1,trial_types_all{:});
end
% if exist('cont_trials_seq', 'var')
%     if sum(max(MMN_trials_seq)> 2)
%         MMN_trials_seq = round(2*MMN_trials_seq/(max(MMN_trials_seq(:))));
%     end
%     
%     MMN_trials_seq_indexed = zeros(size(MMN_trials_seq));
% 
%     stim_params.red_dev_stim_patten = [mode(MMN_trials_seq(:,1)), 3-mode(MMN_trials_seq(:,1));
%                                        mode(MMN_trials_seq(:,2)), 3-mode(MMN_trials_seq(:,2))];
% 
%     MMN_trials_seq_indexed(MMN_trials_seq(:,1) == stim_params.red_dev_stim_patten(1,2),1) = 170;
%     MMN_trials_seq_indexed(MMN_trials_seq(:,2) == stim_params.red_dev_stim_patten(2,2),2) = 270;
% 
%     % mark MMN and flipMMN
%     for n_flip = 1:2
%         temp_seq = 0;
%         for n_trial = 1:size(MMN_trials_seq,1)
%             if  MMN_trials_seq(n_trial,n_flip) == stim_params.red_dev_stim_patten(n_flip,1)
%                 MMN_trials_seq_indexed(n_trial,n_flip) = 1 + n_flip*100 + temp_seq;
%                 temp_seq = temp_seq + 1;
%             else
%                 temp_seq = 0;
%             end
%         end
%     end
% 
%     data.cont_trials_seq = cont_trials_seq;
%     data.MMN_trials_seq = MMN_trials_seq_indexed;
%     data.trial_types = [cont_trials_seq; MMN_trials_seq_indexed(:,1); MMN_trials_seq_indexed(:,2)];

end