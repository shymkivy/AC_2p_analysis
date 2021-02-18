clear;
close all;

addpath('C:\Users\ys2605\Desktop\stuff\AC_2p_analysis\general_functions');

path_csv = 'C:\Users\ys2605\Desktop\stuff\AC_data\behavior\volt\nm_day18_rl_ammn-001\';
fname_csv = 'nm_day18_rl_ammn-001_Cycle00001_VoltageRecording_001.csv';

path_mat = 'C:\Users\ys2605\Desktop\stuff\AC_data\behavior\stim_out\';
fname_mat = 'nm_day18_RL_ammn_1_2_18_21_12h_33m.mat';

data = readtable([path_csv fname_csv]);
data_mat = load([path_mat fname_mat]);

[onset, offset] = f_get_pulse_times(data.x4_LED,.2, 5);

num_trials = data_mat.trial_data.num_trials/2;
num_rewards = data_mat.trial_data.num_rewards;

shift1 = onset(1) - 1;
scale1 = (onset(2) - onset(1))/(data_mat.trial_data.time_paradigm_end);

time_trial_start = data_mat.trial_data.time_trial_start;
time_trial_start(time_trial_start == 0) = [];
time_trial_start = round(time_trial_start*scale1 + shift1);

time_reward_period_start = data_mat.trial_data.time_reward_period_start;
time_reward_period_start(time_reward_period_start == 0) = [];
time_reward_period_start = round(time_reward_period_start*scale1 + shift1);

time_correct_lick = data_mat.trial_data.time_correct_lick;
time_correct_lick(time_correct_lick == 0) = [];
time_correct_lick = round(time_correct_lick*scale1 + shift1);

time_paradigm_end = round(data_mat.trial_data.time_paradigm_end*scale1 + shift1);


data_lick = medfilt1(data.x2_Lick, 5);

data_lick_n = round(data_lick/max(data_lick));
data_lick_on = [0; diff(data_lick_n)]>0.5;

figure; hold on;
plot(data_lick_n)
plot(data_lick_on)
%%


%figure; hold on;

time_onset = zeros(numel(time_trial_start),1);
lick_rate = zeros(numel(time_trial_start),1);

for n_tr = 1:numel(time_trial_start)
    trial_period = (time_trial_start(n_tr)):(time_reward_period_start(n_tr)+5000);
    
    num_licks = sum(data_lick_on(trial_period));
    num_ms = time_reward_period_start(n_tr) - time_trial_start(n_tr);
    
    time_onset(n_tr) = num_ms/1000;
    lick_rate(n_tr) = num_licks/num_ms*1000;
    
    %plot(data_lick(trial_period))
    %plot(data.x3_Stim_type(trial_period))
    %plot(data_lick_on(trial_period))
    %plot(cumsum(data_lick_on(trial_period))./x1', 'color', [.5 .5 .5])
end

time_onset_r = round(time_onset);

time_onset_run = unique(time_onset_r);

thresh = .33;
thresh1 = zeros(numel(time_onset_run),1);
for ii = 1:numel(time_onset_run)
    temp_data = lick_rate(time_onset_r == time_onset_run(ii));
    thresh1(ii) = prctile(temp_data, thresh*100);
end

figure; hold on;
plot(time_onset_r, lick_rate, 'o');
plot(time_onset_run, thresh1);


%%

figure; hold on;
plot(data.Time_ms_, data_lick);
plot(data.Time_ms_, data.x6_LEDBhv);
plot(data.Time_ms_, data.x7_Water);
plot(data.Time_ms_, data.x3_Stim_type);
plot(data.Time_ms_, data.x4_LED);
plot(time_trial_start, 5*ones(numel(time_trial_start),1), '*');
plot(time_reward_period_start, 5*ones(numel(time_reward_period_start),1), '*');
plot(time_correct_lick, 5*ones(numel(time_correct_lick),1), 'k*');
plot(time_paradigm_end, 3*ones(numel(time_paradigm_end),1), '*');
legend('lick', 'LED', 'Water','Stim type', 'synch pulse', 'trial start', 'reward period start', 'correct lick', 'paradigm end')


reward_start_resp = squeeze(f_get_stim_trig_resp(data_lick', time_reward_period_start, [5000 5000]))';
figure; plot(mean(reward_start_resp))
figure; imagesc(reward_start_resp)

trial_start_resp = squeeze(f_get_stim_trig_resp(data_lick', time_trial_start, [3000 5000]))';
figure; plot(mean(trial_start_resp))
figure; imagesc(trial_start_resp)

