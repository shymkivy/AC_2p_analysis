function stats_out = f_dv_stats_within_core(app, params)

peak_stats = params.stats.stat_method; % 'shuff_pool', 'shuff_locwise', 'z_thresh'
stat_source = params.stats.stat_source; % 'All', 'Freqs', 'Freqs_dd'
peak_bin_time = params.stats.peak_bin_time; % sec, .250
num_samp = params.stats.num_shuff_samp; % 1000 before
stat_window = [-2, 6];%params.stats.base_resp_win;
z_thresh = params.stats.z_thresh;
loco_thresh_prc = params.stats.loco_thresh;

n_pl = params.n_pl;
ddata = params.ddata;

fr = 1000/double(ddata.proc_data{1}.frame_data.volume_period);

%%
stat_window_t = (ceil(stat_window(1)*fr):floor(stat_window(2)*fr))/fr;
stat_window_num_baseline_resp_frames = [sum(stat_window_t<=0) sum(stat_window_t>0)];   

%%
win1 = stat_window_num_baseline_resp_frames;
firing_rate = cat(1,params.cdata.S_sm);

[num_cells, T] = size(firing_rate);

%%
num_cells = params.cdata.num_cells;

stim_times_frame_all = ddata.proc_data{1,1}.stim_times_frame;

stim_times = stim_times_frame_all{strcmpi(ddata.proc_ops{1,1}.chan_labels, 'stim type'),n_pl};
stim_times_rew = stim_times_frame_all{strcmpi(ddata.proc_ops{1,1}.chan_labels, 'reward'),n_pl};
stim_times_lick = stim_times_frame_all{strcmpi(ddata.proc_ops{1,1}.chan_labels, 'lick'),n_pl};

num_stim = numel(stim_times);

figure; hold on
plot(stim_times, '.')
plot(stim_times_rew, '.')

stim_is_rev = false(num_stim, 1);
stim_rew_delay = zeros(num_stim,1);
for n_st = 1:num_stim
    temp1 = (stim_times_rew - stim_times(n_st))/fr;
    temp2 = temp1(temp1>0);
    temp3 = temp2 < 3;
    if sum(temp3)
        stim_is_rev(n_st) = 1;
        stim_rew_delay(n_st) = temp2(temp3);
    end
end

%%
trial_types = ddata.trial_types{1};
trial_data_sort = f_get_stim_trig_resp(firing_rate, stim_times, win1);

num_win = sum(win1);
num_trials = numel(stim_times);

%%
samp_all = randsample((win1(1)+1):(T-win1(2)-1), num_samp*num_cells*num_trials, true);
samp_all = reshape(samp_all, num_cells, num_samp * num_trials);
%%
trial_data_samp = zeros(num_cells, num_win, num_samp);
for n_cell = 1:num_cells
    temp_samp = reshape(f_get_stim_trig_resp(firing_rate(n_cell,:), samp_all(n_cell,:), win1), 1, num_win, num_samp, num_trials);
    trial_data_samp(n_cell, :, :) = mean(temp_samp,4);
end

%%
trial_data_samp_mean = mean(trial_data_samp,3);
trial_data_samp_std = std(trial_data_samp, [] ,3);

trial_data_mean = mean(trial_data_sort,3);
trial_data_mean_rew = mean(trial_data_sort(:,:,stim_is_rev),3);
trial_data_mean_notrew = mean(trial_data_sort(:,:,~stim_is_rev),3);

trial_data_mean_z = (trial_data_mean - trial_data_samp_mean)./mean(trial_data_samp_std,2);
trial_data_mean_z_rew = (trial_data_mean_rew - trial_data_samp_mean)./mean(trial_data_samp_std,2);
trial_data_mean_z_notrew = (trial_data_mean_notrew - trial_data_samp_mean)./mean(trial_data_samp_std,2);


stats_out.stim_trial_data_mean = trial_data_mean;
stats_out.trial_data_mean_rew = trial_data_mean_rew;
stats_out.trial_data_mean_notrew = trial_data_mean_notrew;

stats_out.stim_trial_data_mean_z = trial_data_mean_z;
stats_out.trial_data_mean_z_rew = trial_data_mean_z_rew;
stats_out.trial_data_mean_z_notrew = trial_data_mean_z_notrew;

stats_out.num_samp = num_samp;
stats_out.stat_window_t = stat_window_t;
stats_out.stat_window = stat_window;
stats_out.stim_times = stim_times;
stats_out.stim_times_rew = stim_times_rew;
stats_out.stim_times_lick = stim_times_lick;
stats_out.stim_is_rev = stim_is_rev;
stats_out.stim_rew_delay = stim_rew_delay;



%%

[peak_mag, peak_loc] = max(trial_data_mean_z,[],2);
[~, trace_ord] = sort(peak_mag, 'descend');
figure; imagesc(stat_window_t, [], trial_data_mean_z(trace_ord,:))


stim_time_trace = zeros(T,1);
stim_time_rew_trace = zeros(T,1);
stim_time_lick_trace = zeros(T,1);
stim_time_trace(stim_times) = 2;
stim_time_rew_trace(stim_times_rew) = 3;
stim_time_lick_trace(stim_times_lick) = 1;

figure; hold on;
plot(stim_time_trace)
plot(stim_time_rew_trace)
plot(stim_time_lick_trace)


%%
% 
% figure; imagesc(mean(trial_data_samp,3))
% 
% for n_cell_idx = 1:10
%     n_cell = trace_ord(n_cell_idx);
%     trace1_samp = trial_data_samp(n_cell,:,:);
%     mean1_samp = mean(trace1_samp,3);
%     std1_samp = std(trace1_samp, [], 3);
%     trace1 = trial_data_sort(n_cell,:,:);
%     mean1 = mean(trace1,3);
% 
%     figure; hold on;
%     %plot(squeeze(trace1), 'color', [.6 .6 .6]);
%     plot(stat_window_t, mean1_samp, 'k');
%     plot(stat_window_t, mean1_samp + std1_samp*2, 'color', [.6 .6 .6]);
%     plot(stat_window_t, mean1_samp - std1_samp*2, 'color', [.6 .6 .6]);
%     plot(stat_window_t, mean1, 'm');
%     title(sprintf('cell %d', n_cell))
% end

end