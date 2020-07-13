function data = f_s1_process_stim_volt(data, ops)

if ~isfield(data, 'stim_times_trace')
    
    stim_times_trace = cell(1,ops.num_planes);
    for n_pl = 1:ops.num_planes
        % Extract stim onset time
        stim_cutoff = 1/2;
        stim_times_trace{n_pl} = f_s1_get_stim_onsets(data.indexed_volt_data{n_pl}.*logical(data.mmn_phase_mpl{n_pl}),stim_cutoff);
        
        fprintf('Extracted %d stimuli from voltage file for plane %d\n', sum(stim_times_trace{n_pl}), n_pl);
%         figure; hold on;
%         plot(indexed_volt_data{n_pl});
%         plot(stim_times_trace{n_pl});
    end
    data.stim_times_trace = stim_times_trace;
    
end

%%
data.stim_frame_index = cell(1,ops.num_planes);
data.num_trials_extracted = zeros(1,ops.num_planes);
for n_pl = 1:ops.num_planes
    data.stim_frame_index{n_pl} = find(data.stim_times_trace{n_pl});
    data.num_trials_extracted(n_pl) = sum(data.stim_times_trace{n_pl});
end

%%
if ~isfield(data, 'cont_trials_seq')
    % normalize the trigger voltage to index
    num_cont = sum((data.mmn_phase_mpl{1} == 1) .* data.stim_times_trace{1});
    num_mmn = sum((data.mmn_phase_mpl{1} == 2) .* data.stim_times_trace{1});
    
    cont_trials_seq = zeros(num_cont, 1);
    MMN_trials_seq = zeros(num_mmn, 2);
    
    stim_frames = floor(data.stim_params.stim_duration/data.frame_data.volume_period_ave*1000);
    
    % Extract stimulus direction
    %trial_type_directions = zeros(numel(stim_frame_index),1);
    for n_trial = 1:numel(data.stim_frame_index{1})
        % do a median over frames to prevent some unlikely noise
        temp_frame = data.stim_frame_index{1}(n_trial);
        if data.mmn_phase_mpl{1}(data.stim_frame_index{1}(n_trial)) == 1
            cont_trials_seq(n_trial) = median(data.indexed_volt_data{1}(temp_frame:temp_frame+stim_frames));
        elseif data.mmn_phase_mpl{1}(data.stim_frame_index{1}(n_trial)) > 1
            MMN_trials_seq(n_trial - num_cont) = median(data.indexed_volt_data{1}(temp_frame:temp_frame+stim_frames));
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