

path_csv = 'C:\Users\ys2605\Desktop\stuff\AC_data\behavior\volt\nm_day4_RL_ammn_part2-003\';
fname_csv = 'nm_RLammn_day4_part2-003_Cycle00001_VoltageRecording_001.csv';

path_mat = 'C:\Users\ys2605\Desktop\stuff\AC_data\behavior\stim_out\';
fname_mat = 'nm_day4_ready_lick_ammn_part2_2_5_21_14h_46m.mat';

data = readtable([path_csv fname_csv]);

data_mat = load([path_mat fname_mat]);

time_adj_factor = 24*60*60/1e5;

[onset, offset] = f_get_pulse_times(data.x4_LED,.2, 5)

time_trial_start = data_mat.trial_data.time_trial_start*time_adj_factor;
time_reward_period_start = data_mat.trial_data.time_reward_period_start*time_adj_factor;
time_correct_lick = data_mat.trial_data.time_correct_lick*time_adj_factor;
time_paradigm_end = data_mat.trial_data.time_paradigm_end*time_adj_factor;
num_trials = data_mat.trial_data.num_trials/2;
num_rewards = data_mat.trial_data.num_rewards;

time_correct_lick(time_correct_lick==0) = [];
time_correct_lick = round(time_correct_lick*1000);
time_correct_lick_trace = zeros(numel(data.Time_ms_),1);
time_correct_lick_trace(time_correct_lick) = 5;

data_lick = medfilt1(data.x2_Lick, 5);


figure; hold on;
plot(data.Time_ms_, data.x2_Lick);
plot(data.Time_ms_, data_lick);

figure; hold on;
plot(data.Time_ms_, data_lick);
plot(data.Time_ms_, data.x6_LEDBhv);
plot(data.Time_ms_, data.x7_Water);
plot(data.Time_ms_, data.x3_Stim_type);
plot(data.Time_ms_, data.x4_LED);
legend('lick', 'LED', 'Water','Stim type', 'synch pulse')
