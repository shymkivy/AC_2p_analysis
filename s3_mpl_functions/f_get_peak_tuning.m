function peak_tuning_out = f_get_peak_tuning(trial_data_sort, trial_types, trials_to_analyze, trial_window_t, window_time, ops)

[num_cells, ~, num_trials] = size(trial_data_sort);

fr_peak_mag = zeros(num_cells, num_trials);
fr_peak_latency = zeros(num_cells, num_trials);

% extract data from traces
for n_cell = 1:num_cells
    for n_tr = 1:num_trials
        [fr_peak_mag(n_cell,n_tr), fr_peak_latency(n_cell,n_tr)] = max(trial_data_sort(n_cell,:,n_tr));
    end
end

fr_peak_latency_sec = trial_window_t(fr_peak_latency);
fr_peak_latency_sec2 = fr_peak_latency_sec;

fr_peak_latency_sec(fr_peak_latency_sec2<window_time(1)) = NaN;
fr_peak_latency_sec(fr_peak_latency_sec2>window_time(2)) = NaN;

fr_peak_mag(fr_peak_latency_sec2<window_time(1)) = 0;
fr_peak_mag(fr_peak_latency_sec2>window_time(2)) = 0;

%figure; histogram(fr_peak_latency_sec(logical((fr_peak_latency_sec2>window_time(1)).*(fr_peak_latency_sec2<window_time(2)))))

fr_peak_mag_ave = zeros(num_cells, numel(trials_to_analyze));
fr_peak_latency_sec_ave = zeros(num_cells, numel(trials_to_analyze));
for n_tr_type = 1:numel(trials_to_analyze)
    n_tr = trials_to_analyze(n_tr_type);
    fr_peak_mag_ave(:, n_tr_type) = mean(fr_peak_mag(:,trial_types == n_tr),2);
    fr_peak_latency_sec_ave(:, n_tr_type) = nanmean(fr_peak_latency_sec(:,trial_types == n_tr),2);
end   

stat_pk = f_get_stat_thresholds(fr_peak_mag, trial_types, ops.stat.trials_to_sample, trials_to_analyze, ops);

% remove cells that are nonresponsive
remove_cells = stat_pk.sig_thresh(:,1) == 0;

if ~ops.stat.z_scores_average_thresh
    fr_peak_mag_ave_z = (fr_peak_mag_ave - stat_pk.means)./stat_pk.z_factors;
else
    fr_peak_mag_ave_z = (fr_peak_mag_ave - stat_pk.means)./nanmean(stat_pk.z_factors,2);
end
fr_peak_mag_ave_z(remove_cells,:) = min(fr_peak_mag_ave_z(:));
fr_peak_mag_ave(remove_cells,:) = min(fr_peak_mag_ave(:));

%% compute reliability
fr_peak_reliability = zeros(num_cells, n_tr_type);
for n_tr = 1:n_tr_type
    for n_cell = 1:num_cells
        temp_trials = fr_peak_mag(n_cell,trial_types == trials_to_analyze(n_tr));
        fr_peak_reliability(n_cell, n_tr) = sum(temp_trials>(stat_pk.means(n_cell, n_tr)+ops.stat.z_scores_thresh*stat_pk.z_factors(n_cell,n_tr)))/numel(temp_trials);
    end
end

%%
fr_peak_mag_tuned_trials = fr_peak_mag_ave_z>=ops.stat.z_scores_thresh;

%% remove trials
fr_peak_mag_tuned_trials = fr_peak_mag_tuned_trials.*(fr_peak_reliability>=ops.stat.reliability_thresh);

%%
%fr_peak_mag_tuned_trials = fr_peak_mag_ave>stat_pk.sig_thresh;
% fr_peak_mag_tuned_cells_freq = logical(sum(fr_peak_mag_tuned_trials(:,1:10),2));
% [fr_peak_resp_cells_freq_sort_mag, fr_peak_resp_cells_freq_sort_ind] = sort(max(fr_peak_mag_ave_z(:,1:10),[],2), 'descend');
% 
% 
% ctx_tuned_num_cells = zeros(numel(dset_params.ctx_mmn),1);
% ctx_tuned_cells_sort_mag = zeros(num_cells, numel(dset_params.ctx_mmn));
% ctx_tuned_cells_sort_ind = zeros(num_cells, numel(dset_params.ctx_mmn));
% for n_ctx = 1:numel(dset_params.ctx_mmn)
%     ctx_tuned_num_cells(n_ctx) = sum(fr_peak_mag_tuned_trials(:,dset_params.ctx_mmn(n_ctx)));
%     [ctx_tuned_cells_sort_mag(:,n_ctx), ctx_tuned_cells_sort_ind(:,n_ctx)] = sort(fr_peak_mag_ave_z(:,dset_params.ctx_mmn(n_ctx)), 'descend');
% end


peak_tuning_out.fr_peak_mag = fr_peak_mag;
peak_tuning_out.fr_peak_mag_ave = fr_peak_mag_ave;
peak_tuning_out.fr_peak_mag_ave_z = fr_peak_mag_ave_z;
peak_tuning_out.fr_peak_mag_tuned_trials = fr_peak_mag_tuned_trials;
% peak_tuning_out.fr_peak_mag_tuned_cells_freq = fr_peak_mag_tuned_cells_freq;
peak_tuning_out.fr_peak_reliability = fr_peak_reliability;
% peak_tuning_out.fr_peak_resp_cells_freq_sort_mag = fr_peak_resp_cells_freq_sort_mag;
% peak_tuning_out.fr_peak_resp_cells_freq_sort_ind = fr_peak_resp_cells_freq_sort_ind;
% peak_tuning_out.num_freq_tuned_cells = sum(fr_peak_mag_tuned_cells_freq);
% peak_tuning_out.ctx_tuned_num_cells = ctx_tuned_num_cells;
% peak_tuning_out.ctx_tuned_cells_sort_mag = ctx_tuned_cells_sort_mag;
% peak_tuning_out.ctx_tuned_cells_sort_ind = ctx_tuned_cells_sort_ind;
peak_tuning_out.fr_peak_latency_sec = fr_peak_latency_sec;
peak_tuning_out.fr_peak_latency_ave = fr_peak_latency_sec_ave;
peak_tuning_out.stat_pk = stat_pk;

end