function trace_tuning_out = f_get_trace_tuning(trial_data_sort,trial_types, trials_to_analyze,dset_params, ops)

num_cells = size(trial_data_sort,1);

onset_window_frames = dset_params.trial_window{1}.onset_window_frames;
offset_window_frames = dset_params.trial_window{1}.offset_window_frames;
trial_ave_data = f_mpl_trial_average(trial_data_sort,trial_types, trials_to_analyze, 'none');

stat_trace = f_get_stat_thresholds(trial_data_sort, trial_types, ops.stat.trials_to_sample, trials_to_analyze, ops);

tunning_z_mag = (trial_ave_data - stat_trace.means)./nanmean(stat_trace.z_factors,3);

onset_tunning_z_mag = squeeze(max(tunning_z_mag(:,onset_window_frames,:),[],2));
onset_tunning_z_mag(isnan(onset_tunning_z_mag(:))) = 0;
[onset_tunning_sort_z_mag, onset_tunning_sort_ind] = sort(max(onset_tunning_z_mag,[],2), 'descend');

offset_tunning_z_mag = squeeze(max(tunning_z_mag(:,offset_window_frames,:),[],2));
offset_tunning_z_mag(isnan(offset_tunning_z_mag(:))) = 0;
[offset_tunning_sort_z_mag, offset_tunning_sort_ind] = sort(max(offset_tunning_z_mag,[],2), 'descend');

onset_tuned_trials = zeros(num_cells, numel(trials_to_analyze));
onset_tunning_metric = zeros(num_cells, numel(trials_to_analyze));
onset_sensitivity = zeros(num_cells, numel(trials_to_analyze));
for n_cell = 1:num_cells
    for n_trial = 1:numel(trials_to_analyze)
        temp_trace = trial_ave_data(n_cell,onset_window_frames,n_trial);
        temp_thresh = stat_trace.sig_thresh(n_cell,onset_window_frames,n_trial);
        onset_tuned_trials(n_cell,n_trial) = sum(temp_trace(temp_trace == max(temp_trace))>temp_thresh(temp_trace == max(temp_trace)));
        base_mean = mean(mean(trial_ave_data(n_cell,onset_window_frames,(1:numel(trials_to_analyze))~=n_trial),3),2);
        %onset_tunning_mag(n_cell,n_trial) = (mean(temp_trace) - mean(mean(trial_ave_freq(n_cell,onset_window_frames,:),3),2))/(mean(temp_trace) + mean(mean(trial_ave_freq(n_cell,onset_window_frames,:),3),2));
        %onset_tunning_mag(n_cell,n_trial) = (mean(temp_trace))/mean(mean(trial_ave_freq(n_cell,onset_window_frames,:),3),2);
        onset_tunning_metric(n_cell,n_trial) = (mean(temp_trace)-base_mean)/(mean(temp_trace)+base_mean);
        temp_data = squeeze(trial_data_sort(n_cell,onset_window_frames, trial_types == trials_to_analyze(n_trial)));
        onset_sensitivity(n_cell, n_trial) = sum(logical(sum(temp_data>(temp_thresh'))))/sum(trial_types == trials_to_analyze(n_trial));
    end
end
onset_tunning_metric(isnan(onset_tunning_metric)) = 0;
onset_tuned_cells = logical(sum(onset_tuned_trials,2));


offset_tuned_trials = zeros(num_cells, numel(trials_to_analyze));
offset_tunning_metric = zeros(num_cells, numel(trials_to_analyze));
offset_sensitivity = zeros(num_cells, numel(trials_to_analyze));
for n_cell = 1:num_cells
    for n_trial = 1:numel(trials_to_analyze)
        temp_trace = trial_ave_data(n_cell,offset_window_frames,n_trial);
        temp_thresh = stat_trace.sig_thresh(n_cell,offset_window_frames,n_trial);
        offset_tuned_trials(n_cell,n_trial) = sum(temp_trace(temp_trace == max(temp_trace))>temp_thresh(temp_trace == max(temp_trace)));
        base_mean = mean(mean(trial_ave_data(n_cell,offset_window_frames,(1:numel(trials_to_analyze))~=n_trial),3),2);
        offset_tunning_metric(n_cell,n_trial) = (mean(temp_trace)-base_mean)/(mean(temp_trace)+base_mean);
        temp_data = squeeze(trial_data_sort(n_cell,offset_window_frames, trial_types == trials_to_analyze(n_trial)));
        offset_sensitivity(n_cell, n_trial) = sum(logical(sum(temp_data>(temp_thresh'))))/sum(trial_types == trials_to_analyze(n_trial));
    end
end
offset_tunning_metric(isnan(offset_tunning_metric)) = 0;
offset_tuned_cells = logical(sum(offset_tuned_trials,2));

trace_tuning_out.trial_ave = trial_ave_data;
trace_tuning_out.stat_trace = stat_trace;
trace_tuning_out.onset_tuned_trials = onset_tuned_trials;
trace_tuning_out.onset_tuned_cells = onset_tuned_cells;
trace_tuning_out.num_onset_tuned_cells = sum(onset_tuned_cells);
trace_tuning_out.onset_tunning_metric = onset_tunning_metric;
trace_tuning_out.onset_sensitivity = onset_sensitivity;
trace_tuning_out.onset_tunning_z_mag = onset_tunning_z_mag;
trace_tuning_out.onset_tunning_sort_z_mag = onset_tunning_sort_z_mag;
trace_tuning_out.onset_tunning_sort_ind = onset_tunning_sort_ind;
trace_tuning_out.offset_tuned_trials = offset_tuned_trials;
trace_tuning_out.offset_tuned_cells = offset_tuned_cells;
trace_tuning_out.num_offset_tuned_cells = sum(offset_tuned_cells);
trace_tuning_out.offset_tunning_metric = offset_tunning_metric;
trace_tuning_out.offset_sensitivity = offset_sensitivity;
trace_tuning_out.offset_tunning_z_mag = offset_tunning_z_mag;
trace_tuning_out.offset_tunning_sort_z_mag = offset_tunning_sort_z_mag;
trace_tuning_out.offset_tunning_sort_ind = offset_tunning_sort_ind;

end