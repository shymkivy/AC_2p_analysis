function resp = f_get_stim_trig_resp(data, stim_times, trig_window)
% inputs:       
%   data(cell x time_bins)              - traces of signal input
%   stim_times                          - numbers of frames
%   trig_window [baseline, stim_time]   - window around the trigger to extract in frames
%
% outputs
%   resp(cell x resp_window x trials)

[num_cells, ~] = size(data);
num_trials = numel(stim_times);
win_size = sum(trig_window);

resp = zeros(num_cells, win_size, num_trials);

for n_trial = 1:num_trials
    cur_frame = round(stim_times(n_trial));
    time_trace = (cur_frame-trig_window(1)):(cur_frame+trig_window(2)-1);
    resp(:,:,n_trial) = data(:,time_trace);    
end

end